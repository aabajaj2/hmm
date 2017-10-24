import sys
import re
import numpy as np

def take_log(pm):
    '''
    Takes the natural log of Probability Matrices
    :param pm: Probability Matrix
    '''
    for i, row in enumerate(pm):
        for j, item1 in enumerate(row):
            pm[i][j] = np.log(pm[i][j])

def column_epm(nu):
    '''
    gives the column number from epm matrix according to a nucleotide occurrence
    :param nu: single nucleotide
    :return: column number for the emission probability matrix
    '''
    column_in_epm = 0
    if nu == 'a' or nu == 'A':
        column_in_epm = 0
    elif nu == 'c' or nu == 'C':
        column_in_epm = 1
    elif nu == 't' or nu == 'T':
        column_in_epm = 2
    elif nu == 'g' or nu == 'G':
        column_in_epm = 3
    return column_in_epm

def initialize(alpha_matrix, nucleotides, pi_vector, epm, n):
    '''
    Initializes the alpha matrix with pi vector and emission probability matrix
    :param alpha_matrix: The viterbi matrix to be made
    :param nucleotides: the columnns for the gene to analyze
    :param pi_vector: Pi vector for initialization
    :param epm: Emission Probability Matrix
    :param n: number of states
    '''
    one = column_epm(nucleotides[0])
    for i in range(0, n):
        alpha_matrix[i][0] = pi_vector[i] + epm[i][one]

def return_best_path(stack):
    '''
    Converts the best path array to named states, eg. S1, S2, S3
    :param stack: the best path stack
    :return: best path in words
    '''
    best_path = []
    for value in stack:
        state = "S"+str(value+1)
        best_path.append(state)
    return best_path

def read_file():
    '''
    Reads the FASTA file and returns a single string
    :return: nucleotides from the FASTA file
    '''
    filename = sys.argv[1]
    file = open(filename, "r")
    readfile = file.read()
    nucleotides = re.sub('[^a-zA-Z]+', '', readfile)
    return nucleotides

def viterbi_matrix_with_backtrack(nucleotides, alpha_matrix, epm, tpm, n):
    """
    Creates an alpha matrix using viterbi algorithm and also does backtracking to determine the best path with
    Maximum probability
    :param nucleotides: the string of nucleotides returned from the file
    :param alpha_matrix: the result matrix
    :param epm: emission probability matrix
    :param tpm: transition probability matrix
    :param n: number of states
    :return: the best path after back tracking
    """
    backtrack = []
    for iter1 in range(1, len(nucleotides)):
        foralgo = nucleotides[iter1]
        col = column_epm(foralgo)
        sublist = []
        for row in range(n):
            formax = []
            for mx in range(0, n):
                formax.append(alpha_matrix[mx][iter1 - 1] + tpm[mx][row])
            x = max(formax)
            sublist.append(formax.index(x))
            alpha_matrix[row][iter1] = epm[row][col] + x
        backtrack.append(sublist)
    max_probability = alpha_matrix[row][iter1]
    print("*****Maximum Probability*****", max_probability)
    sequence = []
    colnum = len(alpha_matrix[0])-1  #last column values
    for row in range(n):
        sequence.append(alpha_matrix[row][len(alpha_matrix[0])-1])
    index = sequence.index(max(sequence))
    stack = [index]
    while colnum > 0:
        stack.insert(0, backtrack[colnum-1][index])
        index = backtrack[colnum-2][index]
        colnum = colnum - 1
    bestPath = return_best_path(stack)
    return bestPath

def main():
    '''
    Main class for this program
    :return: best path for HMM using viterbi best path algo
    '''
    n = 3
    m = 4
    # can be made command line
    pi_vector = [np.log(0.3), np.log(0.2), np.log(0.5)]
    tpm = [[0.5, 0, 0.5], [0.25, 0.5, 0.25], [0.1, 0.4, 0.5]]
    epm = [[0.3, 0.1, 0.4, 0.2], [0.1, 0.5, 0.1, 0.3], [0.25, 0.25, 0.25, 0.25]]

    # From the example in text
    # pi_vector = [np.log2(0.25), np.log2(0.5), np.log2(0.25)]
    # tpm = [[0.25, 0.5, 0.25], [0.25, 0.25, 0.5], [0.5, 0.5, 0]]
    # epm = [[1, 0, 0, 0], [0.25, 0.5, 0, 0.25], [0.25, 0.25, 0.25, 0.25]]

    take_log(tpm)
    take_log(epm)
    nucleotides = read_file()
    alpha_matrix = np.zeros((n, len(nucleotides)))
    initialize(alpha_matrix, nucleotides, pi_vector, epm, n)
    print("*****Pi Vector*****")
    print(pi_vector)
    print("*****Transition Matrix*****")
    print(tpm)
    print("*****Emission Matrix*****")
    print(epm)
    bestPath = viterbi_matrix_with_backtrack(nucleotides, alpha_matrix, epm, tpm, n)
    print("*****Alpha Matrix:*****")
    print(alpha_matrix)
    print("*******************")
    print("*****Best Path*****")
    print(bestPath)

if __name__ == '__main__':
    main()   