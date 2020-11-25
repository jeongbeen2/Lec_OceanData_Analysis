### SVD script for Week 13 (image process) 
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import os
im_orig = Image.open('marinesnow3.jpg')
im_bw = im_orig.convert('L')
im = np.array(im_bw)
U, s, V = np.linalg.svd(im)

sm = np.zeros([len(U), len(V)])
for j in range(0, len(U)):
    sm[j,j] = s[j]
    
##

small = 40

# U,s,V의 성분을 small만큼 줄이기
Us = U[:,0:small]
sms = sm[0:small, 0:small]
Vs = V[0:small,:]


Us_matmul = np.matmul(Us,sms)

im_recon = np.matmul(Us_matmul, Vs)

plt.imshow(im_recon, cmap='gray')