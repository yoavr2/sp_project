import mykmeanssp as kpp
import numpy as np
import pandas as pd
import sys
import os

import enum


class goal(enum.Enum):
    spk = 1
    wam = 2
    ddg = 3
    lnorm = 4
    jacobi = 5


def calc_DL(row, centroids, cor_num):
    min_diff = sys.maxsize
    for centroid in centroids:
        diff = 0
        for j in range(cor_num):
            diff += ((row[j] - centroid[j]) ** 2)
        if diff < min_diff:
            min_diff = diff
    return min_diff


def kmeans_plus_plus(k, max_iter, eps, file_name):
    with open(file_name) as f:
        content = f.read()
        cor_num = len(content.split("\n")[0].split(","))
        N = len(content.split("\n"))

    lst_centroids = []
    centroids_indices = []
    nodes = pd.read_csv(file_name, header=None, index_col=False)
    np.random.seed(0)
    indexes = [i for i in range(N - 1)]
    index = np.random.choice(indexes)
    centroid = nodes.loc[index]
    lst_centroids.append(centroid)
    centroids_indices.append(index)

    for i in range(1, k):
        nodes['D'] = nodes.apply(lambda row: calc_DL(row, lst_centroids, cor_num), axis=1)
        sum_D = nodes['D'].sum()

        nodes['P'] = nodes['D'].apply(lambda dl: dl / sum_D)

        key = np.random.choice(nodes.index.values, p=nodes['P'])
        lst_centroids.append(nodes.loc[key])
        centroids_indices.append(key)

    print(nodes)
    print(cor_num)
    print(centroids_indices)

    kpp_res = kpp.fit(k, max_iter, eps, nodes.values.tolist(), centroids_indices)   # result will be written to output

    if kpp_res != 0:
        print("An Error Has Occurred")
        sys.exit()

    with open("my_nice_output.txt") as outputf:
        outputf_read = outputf.read()
        the_lines_to_print = outputf_read.split("\n")
        for element in the_lines_to_print:
            print(element)

def cal_jacobi(k, file_name):
    pass


def cal_matrix(k, goal, file_name):
    pass


if __name__ == '__main__':
    """
    args = sys.argv  # args[0] = k, args[1] = goal, args[2] = file_name
    try:
        fr = open("tmpFile.txt", "r")
    except OSError:
        print("Invalid Input!\n")

    try:
        inp_txt = fr.read()
        data_points = inp_txt.split("\n")
        data_num = len(data_points)
        tmp = data_points[0].split(",")
        dimension = len(tmp)
    except:
        print("An Error Has Occurred\n")

    """

    #data_in_float = [[float(x) for x in data_points[i].split(",")] for i in range(data_num-1)]

    data_in_float2 = pd.read_csv("tmpFile.txt", header=None)
    kmeans_plus_plus(2, 300, 0, "tmpFile.txt")

    mat = [[], []]
    mat[0] = [1, 2, 3]
    mat[1] = [4, 5, 6]
    # mat[2] = [1, 1, 1]
    #print(data_in_float)
    #a = kpp.heuristic(data_in_float)
    #print(data_in_float2.values.tolist())
    #a = kpp.jacobi(data_in_float2.values.tolist())
    #print("mat: ")
    #print(a)
    """

    try:
        if len(args) == 3:
            if args[2] == goal.spk:
                kmeans_plus_plus(args[0], 100, 0, args[2])
            if args[2] == goal.jacobi:
                cal_jacobi(args[0], args[2])
            else:
                cal_matrix(args[0], args[1], args[2])
        else:
            print("Invalid Input!\n")
    except:
        pass
        """
