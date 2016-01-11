import scipy.sparse as sps
import numpy as np


def main():
    # index of the row to become the identity row
    rowind = 2

    A = np.arange(16).reshape((4, 4))
    A = sps.csr_matrix(A)
    print 'sparse mat A was before'
    print A.todense()

    print '\nindex of the row that is modified: {0}'.format(rowind)

    setcsrrow2id(A, rowind)
    print '\nsparse mat A is now'
    print A.todense()


def setcsrrow2id(amat, rowind):
    indptr = amat.indptr
    values = amat.data
    indxs = amat.indices

    # get the range of the data that is changed
    rowpa = indptr[rowind]
    rowpb = indptr[rowind+1]

    # new value and its new rowindex
    values[rowpa] = 1.0
    indxs[rowpa] = rowind

    # number of new zero values
    diffvals = rowpb - rowpa - 1

    # filter the data and indices and adjust the range
    values = np.r_[values[:rowpa+1], values[rowpb:]]
    indxs = np.r_[indxs[:rowpa+1], indxs[rowpb:]]
    indptr = np.r_[indptr[:rowind+1], indptr[rowind+1:]-diffvals]

    # hard set the new sparse data
    amat.indptr = indptr
    amat.data = values
    amat.indices = indxs

if __name__ == '__main__':
    main()
