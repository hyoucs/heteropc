import time
import io
from ADMM_TVNL_SAC import ADMM_fully_distributed
import numpy as np
import networkx as nx
import pprint
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import *
import random
from snapvx import *
import scipy.io as sio
import argparse, os
import operator
import json
import matplotlib.pyplot as plt
from numpy import linalg as LA
import cvxpy as cvx
import pickle
import codecs


np.random.seed(1)

Lambda0 = 1.0
Lambda0_1 = 1.0
Lambda1 = 5.0
mu = 0.48
p_norm = 3

rho = 4.0


def main(args):

    inner_iter = args.max_inner_iter

    if args.sac_active:
        iter = args.max_iter
        dim = 4

        graph = nx.read_gpickle(args.graph_dir)
        print(nx.info(graph))
        num_nodes = len(graph.nodes())

        # Bring the features/labels
        with open(args.data_dir, 'r') as fp:
            houses_dict = json.load(fp)
        with open(args.neighbors_dir, 'r') as fp:
            neighbords_dict = json.load(fp)
        with open(args.test_features_dir, 'r') as fp:
            test_features_dict = json.load(fp)

        # Bring the labels
        y = [0] * num_nodes
        W = np.random.rand(num_nodes, (dim - 1))

        for i in range(num_nodes):
            y[i] = float(houses_dict[str(i)][0])
            for k in range(dim - 1):
                W[i][k] = float(houses_dict[str(i)][k + 1])


    if args.cvx_solver_active:

        results = {}

        if args.multiple_runs:
            Lambda1_list = [2.0, 3.0, 4.0] #list(np.arange(2.0, 10.0, 1.0))
            mu_list = [0.44,0.445,0.45,0.455,0.46]#list(np.arange(0.44, 0.5, 0.005))
        else:
            Lambda1_list = [1.0]
            mu_list = [0.5]
        for Lambda1 in Lambda1_list:
            Lambda1 = round(Lambda1,3)
            print('LAMBDA1: ', Lambda1)
            results[Lambda1] = {'nl': {}, 'tvnl': {}}

            nl_solution, variables = cvx_solver_nl(y, W, dim, graph, p_norm, Lambda0, Lambda0_1, Lambda1,
                                                   args.cvx_output_dir, houses_dict)
            print('NL_CVX solver finished with value: ', str(nl_solution.value))
            if args.sac_active:
                nl_total_loss = sac_cvx_test_loss_nl(variables, neighbords_dict, test_features_dict, dim)
                print('NL test loss: ', nl_total_loss)
            else:
                nl_total_loss = cvx_test_loss(variables, neighbords_dict, node_index_dict, test_features_dict, dim)
                print('NL test loss: ', nl_total_loss)

            results[Lambda1]['nl']['train'] = nl_solution.value
            results[Lambda1]['nl']['test'] = nl_total_loss
            if Lambda1 == 0.0:
                continue
            for mu in mu_list:
                mu = round(mu, 3)
                print('MU: ', mu)
                if args.sac_tvnl_active:
                    try:
                        A = ADMM_fully_distributed(y, W, graph, rho, Lambda0, Lambda0_1, Lambda1, mu,
                                                   houses_dict, neighbords_dict, test_features_dict,
                                                   args.tvnl_output_dir, p_norm)
                        tvnl_solution, tvnl_total_loss, tvnl_model = A.runADMM_Grid(iter, inner_iter)
                        print('TVNL_CVX solver finished with value: ', str(tvnl_solution))
                        print('TVNL test loss: ', tvnl_total_loss)
                        results[Lambda1]['tvnl'][mu] = {}
                        results[Lambda1]['tvnl'][mu]['train'] = tvnl_solution
                        results[Lambda1]['tvnl'][mu]['test'] = tvnl_total_loss
                        continue
                    except:
                        continue

                # try:
                #     tvnl_solution, variables = cvx_solver_tvnl(y, W, dim, graph, p_norm, Lambda0, Lambda0_1, Lambda1, mu, args.cvx_output_dir, houses_dict)
                #     print('TVNL_CVX solver finished with value: ', str(tvnl_solution.value))
                #     if args.sac_active:
                #         tvnl_total_loss = sac_cvx_test_loss_tvnl(variables, neighbords_dict, test_features_dict, dim, mu)
                #         print('TVNL test loss: ', tvnl_total_loss)
                #     else:
                #         tvnl_total_loss = cvx_test_loss(variables, neighbords_dict, node_index_dict, test_features_dict, dim)
                #         print('TVNL test loss: ', tvnl_total_loss)
                #
                #     results[Lambda1]['tvnl'][mu] = {}
                #     results[Lambda1]['tvnl'][mu]['train'] = tvnl_solution.value
                #     results[Lambda1]['tvnl'][mu]['test'] = tvnl_total_loss
                # except:
                #     continue
            write_dictionary_to_json(data=results, output_file='Results/result_batch_admm_05.json')

        write_dictionary_to_json(data=results, output_file='Results/result_batch_admm_05.json')


    if args.housing_nl_active:

        gvx = TGraphVX()
        num_nodes = nx.number_of_nodes(graph)
        for i in range(num_nodes):
            x = Variable(dim, name='x')
            obj = 0
            obj += Lambda0 * square(y[i] - ((np.array(W[i, :]) * x[:dim-1]) + x[dim-1]))
            obj += Lambda0_1 * square(norm(x))
            gvx.AddNode(i, Objective=obj)

        for (e0, e1, w) in graph.edges(data=True):
            gvx.AddEdge(e0, e1)
            # gvx.AddEdge(e1, e0)

        gvx.AddEdgeObjectives(laplace_reg)

        start = time.time()
        gvx.Solve(Verbose=True, Rho=rho, MaxIters=100)
        ellapsed = time.time() - start
        print(ellapsed, "seconds; with ADMM")
        # gvx.PrintSolution()

        dif_norm_btw_clusters = 0.0
        num_of_edges_btw_clusters = 0
        dif_norm_within_clusters = 0.0
        num_of_edges_within_clusters = 0

        for (e0, e1) in graph.edges(data=False):
            if houses_dict[index_node_dict[str(e0)]][6] == houses_dict[index_node_dict[str(e1)]][6]:
                dif_norm_within_clusters += LA.norm(np.array(gvx.node_values[e0]) - np.array(gvx.node_values[e1])) / max(
                    LA.norm(gvx.node_values[e0]), LA.norm(gvx.node_values[e1]))
                num_of_edges_within_clusters += 1
            else:
                dif_norm_btw_clusters += LA.norm(np.array(gvx.node_values[e0]) - np.array(gvx.node_values[e1])) / max(
                    LA.norm(gvx.node_values[e0]), LA.norm(gvx.node_values[e1]))
                num_of_edges_btw_clusters += 1

        print('(NL) Within: ', dif_norm_within_clusters / num_of_edges_within_clusters, num_of_edges_within_clusters)
        print('(NL) Between: ', dif_norm_btw_clusters / num_of_edges_btw_clusters, num_of_edges_btw_clusters)
        NL_test_loss = housing_test_loss(gvx, neighbords_dict, node_index_dict, test_features_dict, dim)
        print('(NL) Test Loss: ', NL_test_loss)

        write_variables_into_file(args.nl_output_dir, gvx=gvx, num_nodes=num_nodes, dim=dim)

    if args.karate_nl_active:
        gvx = NL_solver_karate(graph=graph, dim=args.karate_dim, W=W, y=y, write_to_file=True)
        print('Test loss: ', test_Loss(gvx=gvx, W_test=W_test, y_test=y_test, num_nodes=len(graph.nodes()), num_of_observations=args.num_of_observations))

    if args.tvnl_active:
        plotter_y_total, plotter_y_function, plotter_y_reg_func, plotter_y_test, plotter_y_edge , plotter_y_reg_tol = A.runADMM_Grid(iter, inner_iter, p_norm)

        # Prepare plot material
        train_plotter(plotter_y_total=plotter_y_total, plotter_y_function=plotter_y_function,
                      plotter_y_reg_func=plotter_y_reg_func, plotter_y_test=plotter_y_test,
                      plotter_y_edge=plotter_y_edge , plotter_y_reg_tol=plotter_y_reg_tol,
                      iter=iter, inner_iter=inner_iter, output_plot_dir=args.output_plot_dir)

    # Visualize the karate graph
    #vis_tol(graph, A, p_norm)


def sac_cvx_test_loss_tvnl(variables, neighbords_dict, test_features_dict, dim, mu):

    total_test_loss = 0

    # Solve the new node inference problem here!!
    for key, value in neighbords_dict.items():
        key_node_vector = cvx.Variable(dim, name='x_' + str(key))
        a = cvx.Variable(dim, name='a_' + str(key))
        obj = 0
        for neighbor in value:
            variable_arr = np.array(variables[neighbor[0]].value.ravel()).squeeze()
            obj += mu * cvx.norm(key_node_vector - variable_arr + a) + (1 - mu) * cvx.norm(a, p=p_norm)

        problem = cvx.Problem(cvx.Minimize(obj))
        problem.solve()

        y_test = float(test_features_dict[key][0])
        W_test = test_features_dict[key][1:dim]
        key_node_vector = np.array(key_node_vector.value.ravel()).squeeze()
        diff = np.array(y_test - np.matmul(np.array(W_test), np.array(key_node_vector[0:dim - 1])) + np.array(
            key_node_vector[dim - 1]))
        total_test_loss += admm_norm_2_2(diff)

    return float(total_test_loss)/len(neighbords_dict)

def sac_cvx_test_loss_nl(variables, neighbords_dict, test_features_dict, dim):

    total_test_loss = 0

    # Solve the new node inference problem here!!

    for key, value in neighbords_dict.items():
        key_node_vector = cvx.Variable(dim, name='x_' + str(key))
        obj = 0
        for neighbor in value:
            variable_arr = np.array(variables[neighbor[0]].value.ravel()).squeeze()
            obj += cvx.norm(key_node_vector - variable_arr)

        problem = cvx.Problem(cvx.Minimize(obj))
        problem.solve()

        y_test = float(test_features_dict[key][0])
        W_test = test_features_dict[key][1:dim]
        key_node_vector = np.array(key_node_vector.value.ravel()).squeeze()
        diff = np.array(y_test - np.matmul(np.array(W_test), np.array(key_node_vector[0:dim - 1])) + np.array(
            key_node_vector[dim - 1]))
        total_test_loss += admm_norm_2_2(diff)

    return float(total_test_loss)/len(neighbords_dict)

# Modify the test loss so that it will only have (y-y_pred)**2
def cvx_test_loss(variables, neighbords_dict, node_index_dict, test_features_dict, dim):
    total_test_loss = 0
    for key, value in neighbords_dict.items():
        key_node_vector = np.zeros(dim)
        for neighbor in value:
            for d in range(dim):
                variable_arr = np.array(variables[node_index_dict[neighbor[0]]].value.ravel()).squeeze()
                key_node_vector[d] += variable_arr[d]

        key_node_vector = 1.0 * key_node_vector / len(value)

        y_test = float(test_features_dict[key][0])
        W_test = test_features_dict[key][1:dim]

        diff = np.array(y_test - np.matmul(np.array(W_test), np.array(key_node_vector[0:dim - 1])) + np.array(key_node_vector[dim - 1]))
        total_test_loss += admm_norm_2_2(diff)

    return total_test_loss

def test_Loss(gvx, W_test, y_test, num_nodes, num_of_observations):
    test_loss = 0
    for i in range(num_nodes):
        for j in range(num_of_observations):
            # Add function loss
            diff = np.array(y_test[i,j] - np.matmul(np.array(W_test[i,:,j]), np.array(gvx.node_values[i])))
            test_loss += LA.norm(diff)**2

        # Add regularization loss on X
        #test_loss += Lambda0_1 * LA.norm(np.array(gvx.node_values[i]))**2
    return test_loss

def housing_test_loss(gvx, neighbords_dict, node_index_dict, test_features_dict, dim):
    total_test_loss = 0
    for key, value in neighbords_dict.items():
        key_node_vector = np.zeros(dim)
        for neighbor in value:
            for d in range(dim):
                key_node_vector[d] += gvx.node_values[node_index_dict[neighbor[0]]][d]

        key_node_vector = 1.0 * key_node_vector/len(value)

        y_test = float(test_features_dict[key][0])
        W_test = test_features_dict[key][1:dim]

        diff = np.array(y_test - np.matmul(np.array(W_test), np.array(key_node_vector[0:dim-1])) + np.array(key_node_vector[dim-1]))
        total_test_loss += admm_norm_2_2(diff)

    return total_test_loss


def admm_norm_2_2(input_vector):
    output = 0.0
    for i in np.nditer(input_vector):
        output += i**2

    return output

def vis_tol(nxG, A, p_norm):
    p_norm = p_norm
    partition = []
    color_map = []
    font_color_map = []
    officers = []
    mrhis = []
    for i in range(len(nxG.nodes())):
        if nxG.node[i]['club'] == 'Officer':
            officers.append(i)
            partition.append(0)
            color_map.append('Yellow')
            font_color_map.append('Black')
        else:
            mrhis.append(i)
            partition.append(1)
            color_map.append('Purple')
            font_color_map.append('White')

    pos = nx.spring_layout(nxG)  # compute graph layout

    # val_map = {'0': 1.0,
    #            '1': 0.5714285714285714}
    #
    # values = [val_map.get(node) for node in nxG.nodes()]
    plt.figure(figsize=(8, 8))  # image is 8 x 8 inches
    plt.axis('off')
    nx.draw_networkx_nodes(nxG, pos, node_size=300, node_color=color_map)
    #nx.draw_networkx_nodes(nxG, pos, node_size=300, cmap=plt.cm.RdYlGn, node_color=color_map)
    tol_pnorm = []
    for (e0, e1, w) in nxG.edges(data=True):
        val = LA.norm(A.X[e0,:] - A.X[e1,:],p_norm)/(max(LA.norm(A.X[e0,:],p_norm), LA.norm(A.X[e1,:],p_norm)))
        #val = LA.norm((A.X[e0, :] - A.X[e1, :]),p_norm)

        #val = LA.norm((A.A[e0, e1, :]), p_norm)
        tol_pnorm.append(val)
        nxG[e0][e1]['tol'] = str(round(val,2))
    nx.draw_networkx_edges(nxG, pos, edge_color=tol_pnorm, edge_cmap=plt.cm.Reds, alpha=0.3, width=4)
    # nx.draw_networkx_edge_labels(nxG, pos, labels=nx.get_edge_attributes(nxG,'tol'))
    nx.draw_networkx_labels(nxG.subgraph(officers), pos, font_color='k') #labels=nx.get_node_attributes(nxG,'club')
    nx.draw_networkx_labels(nxG.subgraph(mrhis), pos, font_color='w')
    plt.show(nxG)

    # plt.figure(figsize=(8, 8))  # image is 8 x 8 inches
    # plt.axis('off')
    # nx.draw_networkx_nodes(nxG, pos, node_size=200, cmap=plt.cm.RdYlBu, node_color=partition)
    # nlasso_norm = []
    # for (e0, e1, w) in nxG.edges(data=True):
    #     val = LA.norm(A.U[e0,e1,:]-A.U[e1,e0,:]+A.A[e0,e1,:],2)
    #     nlasso_norm.append(val)
    #     nxG[e0][e1]['tol'] = str(round(val,2))
    # nx.draw_networkx_edges(nxG, pos, edge_color=nlasso_norm, edge_cmap=plt.cm.Reds, alpha=0.3, width=4)
    # # nx.draw_networkx_edge_labels(nxG, pos, labels=nx.get_edge_attributes(nxG,'tol'))
    # plt.show(nxG)

def laplace_reg(src, dst, data):
    return (norm(src['x'] - dst['x']), [])

def generate_karate_graph(dim, num_of_observations, scaler_vector, scaler_observations, scaler_label):

    graph = nx.karate_club_graph()  # Bring the graph here
    print(nx.info(graph))
    num_nodes = len(graph.nodes())

    # Construct ground truth vectors for each community
    g_truth_vectors_0 = scaler_vector * np.random.rand(dim)  # For officer
    g_truth_vectors_1 = scaler_vector * np.random.rand(dim)  # For Mr. Hi

    # Generate observations for each node
    W = scaler_observations * np.random.rand(num_nodes, dim, num_of_observations)  # Bring the features

    # Now construct labels using ground truth vectors
    y = np.random.rand(num_nodes, num_of_observations)

    for i in range(num_nodes):
        for j in range(num_of_observations):
            if graph.node[i]['club'] == 'Officer':
                y[i, j] = np.dot(np.array(W[i, :, j]), np.array(g_truth_vectors_0)) + scaler_label
            else:
                y[i, j] = np.dot(np.array(W[i, :, j]), np.array(g_truth_vectors_1)) + scaler_label

    # Generate test observations and labels as well. This time do not add noise to labels!!
    y_test = np.random.rand(num_nodes, num_of_observations)
    W_test = scaler_observations * np.random.rand(num_nodes, dim, num_of_observations)
    for i in range(num_nodes):
        for j in range(num_of_observations):
            if graph.node[i]['club'] == 'Officer':
                y_test[i, j] = np.dot(np.array(W_test[i, :, j]), np.array(g_truth_vectors_0))
            else:
                y_test[i, j] = np.dot(np.array(W_test[i, :, j]), np.array(g_truth_vectors_1))

    return graph, g_truth_vectors_0, g_truth_vectors_1, W, y, W_test, y_test

def write_variables_into_file(output_dir, gvx, num_nodes, dim):

    with open(output_dir, 'w') as file:
        for i in range(num_nodes):
            for d in range(dim):
                file.write(str(gvx.node_values[i][d]))
                if d == dim - 1:
                    break
                file.write(' ')
            file.write('\n')

def write_labels_into_file(output_dir, graph, num_nodes):
    with open(output_dir, 'w') as file:
        for i in range(num_nodes):
            if graph.node[i]['club'] == 'Officer':
                file.write(str(2))
            else:
                file.write(str(1))
            file.write('\n')

def train_plotter(plotter_y_total, plotter_y_function, plotter_y_reg_func, plotter_y_test, plotter_y_edge,
                  plotter_y_reg_tol, iter, inner_iter, output_plot_dir):

    random_x = np.linspace(0, iter, iter)
    # Create traces
    trace1 = go.Scatter(
        x=random_x,
        y=plotter_y_total,
        mode='lines+markers',
        name='total train loss'
    )
    trace2 = go.Scatter(
        x=random_x,
        y=plotter_y_function,
        mode='lines+markers',
        name='function train loss ' + str(Lambda0)
    )
    trace3 = go.Scatter(
        x=random_x,
        y=plotter_y_reg_func,
        mode='lines+markers',
        name='function reg loss ' + str(Lambda0_1)
    )
    trace4 = go.Scatter(
        x=random_x,
        y=plotter_y_test,
        mode='lines+markers',
        name='function test loss ' + str(Lambda0)
    )
    trace5 = go.Scatter(
        x=random_x,
        y=plotter_y_edge,
        mode='lines+markers',
        name='edge loss ' + str(Lambda1)
    )

    trace6 = go.Scatter(
        x=random_x,
        y=plotter_y_reg_tol,
        mode='lines+markers',
        name='tolerance reg loss ' + str(Lambda2)
    )

    data = [trace1, trace2, trace3, trace4, trace5, trace6]

    layout = dict(title='rho = ' + str(rho) + ' inner_iter = ' + str(inner_iter) + ' p_norm = ' + str(p_norm))

    fig = dict(data=data, layout=layout)

    plotly.offline.plot(fig, filename=output_plot_dir)

def NL_solver_karate(graph, dim, W, y, write_to_file=False):
    gvx = TGraphVX()
    num_nodes = nx.number_of_nodes(graph)
    for i in range(num_nodes):
        x = Variable(dim, name='x')
        obj = 0
        for j in range(args.num_of_observations):
            obj += Lambda0 * square(y[i, j] - (np.array(W[i, :, j]) * x))
        obj += Lambda0_1 * square(norm(x))
        gvx.AddNode(i, Objective=obj)

    for (e0, e1, w) in graph.edges(data=True):
        gvx.AddEdge(e0, e1)
        # gvx.AddEdge(e1, e0)

    gvx.AddEdgeObjectives(laplace_reg)

    start = time.time()
    gvx.Solve(Verbose=True, Rho=rho, MaxIters=100)
    ellapsed = time.time() - start
    print(ellapsed, "seconds; with ADMM")
    # gvx.PrintSolution()
    #
    dif_norm_btw_clusters = 0.0
    num_of_edges_btw_clusters = 0
    dif_norm_within_clusters = 0.0
    num_of_edges_within_clusters = 0

    for (e0, e1) in graph.edges(data=False):
        if graph.node[e0]['club'] == graph.node[e1]['club']:
            dif_norm_within_clusters += LA.norm(np.array(gvx.node_values[e0]) - np.array(gvx.node_values[e1])) / max(
                LA.norm(gvx.node_values[e0]), LA.norm(gvx.node_values[e1]))
            num_of_edges_within_clusters += 1
        else:
            dif_norm_btw_clusters += LA.norm(np.array(gvx.node_values[e0]) - np.array(gvx.node_values[e1])) / max(
                LA.norm(gvx.node_values[e0]), LA.norm(gvx.node_values[e1]))
            num_of_edges_btw_clusters += 1

    print('(NL) Within: ', dif_norm_within_clusters / num_of_edges_within_clusters, num_of_edges_within_clusters)
    print('(NL) Between: ', dif_norm_btw_clusters / num_of_edges_btw_clusters, num_of_edges_btw_clusters)


    if write_to_file:
        write_variables_into_file('Viz/Karate_X_NL_2.txt', gvx=gvx, num_nodes=num_nodes, dim=dim)
        #write_labels_into_file('Viz/Karate_X_labels.txt', graph=graph, num_nodes=num_nodes)

    return gvx


# Save tolerance variables too in order to make further analysis
def cvx_solver_tvnl(y, W, dim, graph, p_norm, Lambda0, Lambda0_1, Lambda1, mu, cvx_output_dir, houses_dict):


    # Create two scalar optimization variables.
    num_nodes = nx.number_of_nodes(graph)
    variables = {i: cvx.Variable(dim, name='X_' + str(i)) for i in range(num_nodes)}
    tolerances = {}

    obj = 0
    for i in range(num_nodes):
        obj += Lambda0 * cvx.square(cvx.norm(y[i] - ((np.array(W[i, :]) * variables[i][:dim - 1]) + variables[i][dim - 1])))
        obj += Lambda0_1 * cvx.square(cvx.norm(variables[i]))
        for Id in graph.neighbors(i):
            if Id > i:
                a = cvx.Variable(dim, name='a_' + str(i) + '_' + str(Id))
                tolerances[str(i)+','+str(Id)] = a
                obj += (float(Lambda1)/mu) * (mu * cvx.norm(variables[i] - variables[Id] + tolerances[str(i)+','+str(Id)]) + (1-mu) * cvx.norm(tolerances[str(i)+','+str(Id)], p=p_norm))
                #obj += Lambda2 * cvx.norm(a, p=p_norm)

    print('Solving for tvnl...')
    problem = cvx.Problem(cvx.Minimize(obj))
    problem.solve()

    with open('Outputs/sac_tolerances_' + cvx_output_dir, 'w') as file:
        for key in tolerances:
            variable_arr = np.array(tolerances[key].value.ravel()).squeeze()
            file.write(key + ':')
            for d in range(dim):
                file.write(str(variable_arr[d]))
                if d == dim - 1:
                    break
                file.write(' ')
            file.write('\n')

    with open('Outputs/sac_tvnl_' + cvx_output_dir, 'w') as file:
        for i in range(num_nodes):
            variable_arr = np.array(variables[i].value.ravel()).squeeze()
            for d in range(dim):
                file.write(str(variable_arr[d]))
                if d == dim - 1:
                    break
                file.write(' ')
            file.write('\n')

    with open('Viz/sac_labels.txt', 'w') as file:
        for i in range(num_nodes):
            file.write(str(houses_dict[str(i)][6]))
            file.write('\n')

    cluster_dif_computer(graph, variables, houses_dict)

    return problem, variables


def cvx_solver_nl(y, W, dim, graph, p_norm, Lambda0, Lambda0_1, Lambda1, cvx_output_dir, houses_dict):

    # Create two scalar optimization variables.
    num_nodes = nx.number_of_nodes(graph)
    variables = {i: cvx.Variable(dim, name='X_' + str(i)) for i in range(num_nodes)}
    obj = 0
    for i in range(num_nodes):
        obj += Lambda0 * cvx.square(cvx.norm(y[i] - ((np.array(W[i, :]) * variables[i][:dim - 1]) + variables[i][dim - 1])))
        obj += Lambda0_1 * cvx.square(cvx.norm(variables[i]))
        for Id in graph.neighbors(i):
            if Id > i:
                obj += Lambda1 * cvx.norm(variables[i] - variables[Id])

    print('Solving for nl...')
    problem = cvx.Problem(cvx.Minimize(obj))
    problem.solve()

    with open('Outputs/sac_nl_' + cvx_output_dir, 'w') as file:
        for i in range(num_nodes):
            variable_arr = np.array(variables[i].value.ravel()).squeeze()
            for d in range(dim):
                file.write(str(variable_arr[d]))
                if d == dim - 1:
                    break
                file.write(' ')
            file.write('\n')

    cluster_dif_computer(graph, variables, houses_dict)

    return problem, variables


def cluster_dif_computer(graph, variables, houses_dict):
    dif_norm_btw_clusters = 0.0
    num_of_edges_btw_clusters = 0
    dif_norm_within_clusters = 0.0
    num_of_edges_within_clusters = 0
    for (e0, e1) in graph.edges(data=False):
        x0 = np.array(variables[e0].value.ravel()).squeeze()
        x1 = np.array(variables[e1].value.ravel()).squeeze()
        if houses_dict[str(e0)][6] == houses_dict[str(e1)][6]:
            dif_norm_within_clusters += admm_norm_2(x0 - x1) / max(
                admm_norm_2(x0), admm_norm_2(x1))
            num_of_edges_within_clusters += 1
        else:
            dif_norm_btw_clusters += admm_norm_2(x0 - x1) / max(
                admm_norm_2(x0), admm_norm_2(x1))
            num_of_edges_btw_clusters += 1

    print('Diff Within:  ', dif_norm_within_clusters / num_of_edges_within_clusters,
          num_of_edges_within_clusters)
    print('Diff Between: ', dif_norm_btw_clusters / num_of_edges_btw_clusters, num_of_edges_btw_clusters)

def admm_norm_2_2(input_vector):
    output = 0.0
    for i in np.nditer(input_vector):
        output += i**2

    return output

def admm_norm_2(input_vector):
    output = 0.0
    for i in np.nditer(input_vector):
        output += i**2

    return math.sqrt(output)


def write_dictionary_to_json(data, output_file):
	with codecs.open(output_file, 'w', 'utf8') as out_writer:
		out_writer.write(json.dumps(data, indent=4, separators=(',', ': '), ensure_ascii=False))

if __name__ == "__main__":

    parser = argparse.ArgumentParser("main.py")

    parser.add_argument('--data_dir', type=str,
                        help='data dir which hold network and features')
    parser.add_argument('--graph_dir', type=str,
                        help='directory pointing to networkx graph')
    parser.add_argument('--index_node_dir', type=str,
                        help='dictionary matching index and node keys')
    parser.add_argument('--node_index_dir', type=str,
                        help='dictionary matching node and index keys')
    parser.add_argument('--neighbors_dir', type=str,
                        help='')
    parser.add_argument('--test_features_dir', type=str,
                        help='')

    parser.add_argument('--tvnl_output_dir', type=str,
                        help='Output dir to save learned features for each node')

    parser.add_argument('--nl_output_dir', type=str,
                        help='Output dir to save learned features for each node')

    parser.add_argument('--cvx_output_dir', type=str,
                        help='Output dir to save learned features for each node')

    parser.add_argument('--generate_random', type=bool, default=True,
                        help='Determines whether to generate random graph/data or take a real world example')

    parser.add_argument('--max_iter', type=int, default=100,
                        help='Maximum number of outer iterations until termination')
    parser.add_argument('--max_inner_iter', type=int, default=20,
                        help='Maximum number of inner iterations (U & A updates)')

    parser.add_argument('--karate_dim', type=int, default=5,
                        help='Number of dimensions for node variables')
    parser.add_argument('--num_of_observations', type=int, default=2,
                        help='Number of observations for each node in karate graph')

    parser.add_argument('--scaler_vector', type=float, default=10.0,
                        help='Scaler for randomly generated vectors for karate graph')
    parser.add_argument('--scaler_observations', type=float, default=10.0,
                        help='Scaler for randomly generated observations for karate graph')
    parser.add_argument('--scaler_label', type=float, default=2.0,
                        help='Scaler for pre-computed labels for karate graph')

    parser.add_argument('--max_karate_iter', type=int, default=100,
                        help='Maximum number of outer iterations for karate graph')

    parser.add_argument('--karate_nl_active', default=False, action='store_true',
                        help='')
    parser.add_argument('--housing_nl_active', default=False, action='store_true',
                        help='')

    parser.add_argument('--tvnl_active', default=False, action='store_true',
                        help='')

    parser.add_argument('--cvx_solver_active', default=False, action='store_true',
                        help='')

    parser.add_argument('--sac_active', default=False, action='store_true',
                        help='')
    parser.add_argument('--sac_tvnl_active', default=False, action='store_true',
                        help='')

    parser.add_argument('--output_plot_dir', type=str,
                        help='Output dir to save losses plot')

    parser.add_argument('--multiple_runs', default=False, action='store_true',
                        help='')

    args = parser.parse_args()
    main(args)