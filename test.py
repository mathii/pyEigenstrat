# Tests for the pyEigenstrat class.
# Test data located in ./testdata

from __future__ import division, print_function
import pyEigenstrat, pdb
import numpy as np

test_gt=np.array([
    [0, 2, 0, 0, 1, 0],
    [0, 0, 9, 9, 9, 1],
    [2, 2, 2, 1, 2, 2],
    [0, 0, 0, 0, 0, 0],
    [2, 2, 2, 2, 2, 2],
    [0, 1, 0, 1, 0, 1],
    [9, 0, 0, 1, 0, 9],
    [0, 1, 2, 1, 0, 1]
], dtype='i1')

print("Starting test")

# 1. test data loading
for what in ["unpacked", "packed"]:
    data=pyEigenstrat.load("testdata/"+what) 
    load_gt=data.geno()
    # Test that loading all the data works
    if (load_gt==test_gt).all():
        print("Loaded data matches test data ("+what+")")
    else:
        raise Exception("Loaded data does not match test data ("+what+")")

    # Test that iterating works
    for i,d in enumerate(data):
        if not (d==test_gt[i,:]).all():
            raise Exception("Loaded data does not match test data ("+what+")")
    print("Iterated data matches test data ("+what+")")

    # Test that we can load specific individuals and snps
    data=pyEigenstrat.load("testdata/"+what, pops=["POP1"], inds=["IND5"], snps=["SNP_1", "SNP_4", "SNP_7", "SNP_8"])
    subtest_gt=test_gt[np.array([0,3,6,7]),:][:,np.array([0,1,3,4])]
    load_gt=data.geno()
    # Test that loading all the data works
    if (load_gt==subtest_gt).all():
        print("Loaded sub-data matches test data ("+what+")")
    else:
        raise Exception("Loaded data does not match test data ("+what+")")

    # Test that iterating works
    for i,d in enumerate(data):
        if not (d==subtest_gt[i,:]).all():
            raise Exception("Loaded data does not match test data ("+what+")")
    print("Iterated sub-data matches test data ("+what+")")

# 2 Test that packed and unpacked files look the same
# packedfiles=["testdata/packed", "testdata_ho/hou3"]
# unpackedfiles=["testdata/unpacked", "testdata_ho/hov3"]
packedfiles=["testdata/packed"]
unpackedfiles=["testdata/unpacked"]
for i in range(len(packedfiles)):
    pack_data=pyEigenstrat.load(packedfiles[i]) 
    unpack_data=pyEigenstrat.load(unpackedfiles[i]) 
    if (pack_data.geno()==unpack_data.geno()).all():
        print("Packed and unpacked data matches (%s,%s)"%(packedfiles[i], unpackedfiles[i]))
    else:
         raise Exception("Packed and unpacked data doesn't match (%s,%s)"%(packedfiles[i], unpackedfiles[i]))

