#import numpy as np
#import pandas as pd
import sys
import os
import mykmeanssp as kpp
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
            diff += ((row[j] - centroid['x_' + str(j)]) ** 2)
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
    nodes = pd.read_csv(file_name, header=None, )
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

    kpp_res = kpp.fit(k, max_iter, eps, "my_nice_output.txt")   # result will be written to output

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

    args = sys.argv  # args[0] = k, args[1] = goal, args[2] = file_name
    a = kpp.test("tmpFile.txt")
    print(a)
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
