import numpy as np
import scipy.io
from scipy.stats import stats
import brainiak.funcalign.srm
import sys

inputfile = sys.argv[1]
iteration = int(sys.argv[2])
feature = int(sys.argv[3])

inputdata = scipy.io.loadmat(inputfile)
origdata = inputdata['origdata']
subjects = origdata.shape[2]
# re-arrange dataarray
data2trans = []
for s in range(subjects):
    x1 = np.delete(origdata[:,:,s], np.where(np.isnan(origdata[:,:,s]))[0], axis=0)
    data2trans.append(x1)
# inputdata normalization
#for subject in range(subjects):
#    data2trans[subject] = stats.zscore(data2trans[subject],axis=1,ddof=1)

srm = brainiak.funcalign.srm.SRM(n_iter=iteration, features=feature)
srm.fit(data2trans)
srmdata = srm.transform(data2trans)

# outputdata normalization
#for subject in range(subjects):
#    srmdata[subject] = stats.zscore(srmdata[subject], axis=1, ddof=1)

outputdata = {'srmdata':[]};
outputdata['srmdata'] = srmdata
scipy.io.savemat('srmData.mat', outputdata)
