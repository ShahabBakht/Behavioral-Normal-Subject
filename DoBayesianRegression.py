# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 09:44:57 2016

@author: labuser
"""

import numpy as np
import scipy as sp

from sklearn import linear_model


listOfSubjects = ['ag', 'az' , 'cs', 'gc', 'hr', 'ls', 'mp', 'sb', 'tc', 'vs']


T = np.linspace(5, 20, 20)[:, np.newaxis];

for subjcount in range(0,9):
    subjectIdx = listOfSubjects[subjcount]
    fileName = 'D:\Analysis\Behavioral-Normal-Subject\Final Clean Data\SVC_prctl_' + subjectIdx + '.mat'
    DATAfromMATLAB = sp.io.loadmat(fileName)
    
    data2 = DATAfromMATLAB.get('prctlSVC_2')
    Xdata2 = data2[0,:];
    Xdata2 = np.matrix(Xdata2).reshape(-1,1);
    Ydata2 = data2[1,:];
    Ydata2 = np.matrix(Ydata2).reshape(-1,1);
    
    data6 = DATAfromMATLAB.get('prctlSVC_6')
    Xdata6 = data6[0,:];
    Xdata6 = np.matrix(Xdata6).reshape(-1,1);
    Ydata6 = data6[1,:];
    Ydata6 = np.matrix(Ydata6).reshape(-1,1);

    data20 = DATAfromMATLAB.get('prctlSVC_20')
    Xdata20 = data20[0,:];
    Xdata20 = np.matrix(Xdata20).reshape(-1,1);
    Ydata20 = data20[1,:];
    Ydata20 = np.matrix(Ydata20).reshape(-1,1);

    clf = linear_model.BayesianRidge(compute_score=True,normalize=True)

    y2_ = clf.fit(Xdata2, Ydata2).predict(T)
    c2 = clf.coef_
    a2 = clf.alpha_
    l2 = clf.lambda_
    Score2 = clf.score(Xdata2,Ydata2)

    y6_ = clf.fit(Xdata6, Ydata6).predict(T)
    c6 = clf.coef_
    a6 = clf.alpha_
    l6 = clf.lambda_
    Score6 = clf.score(Xdata6,Ydata6)

    y20_ = clf.fit(Xdata20, Ydata20).predict(T)
    c20 = clf.coef_
    a20 = clf.alpha_
    l20 = clf.lambda_
    Score20 = clf.score(Xdata20,Ydata20)

    Alpha = [a2,a6,a20]
    Lambda = [l2,l6,l20]
    Score = [Score2,Score6,Score20]
    
    resultName = 'D:\Analysis\Behavioral-Normal-Subject\Final Clean Data\Alpha_Lambda_Score_' + subjectIdx + '.mat'
    sp.io.savemat(resultName,{'Alpha':Alpha,'Lambda':Lambda,'Score':Score})
    
    
    
