"""
============================
Nearest Neighbors regression
============================

Demonstrate the resolution of a regression problem
using a k-Nearest Neighbor and the interpolation of the
target using both barycenter and constant weights.

"""
print(__doc__)

# Author: Shahab Bakhtiari <shahab.bakhtiari@mail.mcgill.ca>



###############################################################################
# Generate sample data
import numpy as np
import matplotlib.pyplot as plt
from sklearn import neighbors
from numpy import genfromtxt
from sklearn.cross_validation import KFold
from sklearn import linear_model

#np.random.seed(0)
#X = np.sort(5 * np.random.rand(40, 1), axis=0)
#T = np.linspace(0, 5, 500)[:, np.newaxis]
#y = np.sin(X).ravel()
#
## Add noise to targets
#y[::5] += 1 * (0.5 - np.random.rand(8))

###############################################################################
# Load data
#Xdata1 = genfromtxt('X1.csv', delimiter=',')
#Xd1 = np.matrix(Xdata1).reshape(-1,1)
#Ydata1 = genfromtxt('Y1.csv', delimiter=',')
#Yd1 = np.matrix(Ydata1).reshape(-1,1)

Xdata6 = genfromtxt('X6.csv', delimiter=',')
Xd6 = np.matrix(Xdata6).reshape(-1,1)
Ydata6 = genfromtxt('Y6.csv', delimiter=',')
Yd6 = np.matrix(Ydata6).reshape(-1,1)

Xdata10 = genfromtxt('X10.csv', delimiter=',')
Xd10 = np.matrix(Xdata10).reshape(-1,1)
Ydata10 = genfromtxt('Y10.csv', delimiter=',')
Yd10 = np.matrix(Ydata10).reshape(-1,1)

Xdata20 = genfromtxt('X20.csv', delimiter=',')
Xd20 = np.matrix(Xdata20).reshape(-1,1)
Ydata20 = genfromtxt('Y20.csv', delimiter=',')
Yd20 = np.matrix(Ydata20).reshape(-1,1)


T = np.linspace(5, 20, 20)[:, np.newaxis];

###############################################################################
# Fit regression model
n_neighbors = 10
knn = neighbors.KNeighborsRegressor(n_neighbors, weights='uniform')
y1_ = knn.fit(Xd1, Yd1).predict(T)

y6_ = knn.fit(Xd6, Yd6).predict(T)

y10_ = knn.fit(Xd10, Yd10).predict(T)

y20_ = knn.fit(Xd20, Yd20).predict(T)

fig1 = plt.figure()
plt.scatter(Xd1, Yd1, c='k', label='data')
line_1, = plt.plot(T, y1_, 'k', label='1 degree')
legends = plt.legend(handles=[line_1], loc=2)
ax = plt.gca().add_artist(legends)
plt.axis([9,21,-5,20])
plt.title("KNeighborsRegressor (k = %i)" % (n_neighbors))
plt.xlabel('Target Velocity (degree/s)')
plt.ylabel('Eye Velocity (degree/s)')
plt.show()

fig2 = plt.figure()
plt.scatter(Xd6, Yd6, c='b', label='data')
line_6, = plt.plot(T, y6_, 'b', label='6 degrees')
legends = plt.legend(handles=[line_6], loc=2)
ax = plt.gca().add_artist(legends)
plt.axis([9,21,-5,20])
plt.title("KNeighborsRegressor (k = %i)" % (n_neighbors))
plt.xlabel('Target Velocity (degree/s)')
plt.ylabel('Eye Velocity (degree/s)')
plt.show()

fig3 = plt.figure()
plt.scatter(Xd10, Yd10, c='g', label='data')
line_10, = plt.plot(T, y10_, 'g', label='10 degrees')
legends = plt.legend(handles=[line_10], loc=2)
ax = plt.gca().add_artist(legends)
plt.axis([9,21,-5,20])
plt.title("KNeighborsRegressor (k = %i)" % (n_neighbors))
plt.xlabel('Target Velocity (degree/s)')
plt.ylabel('Eye Velocity (degree/s)')
plt.show()

fig4 = plt.figure()
plt.scatter(Xd20, Yd20, c='r', label='data')
line_20, = plt.plot(T, y20_, 'r', label='20 degrees')
legends = plt.legend(handles=[line_20], loc=2)
ax = plt.gca().add_artist(legends)
plt.axis([9,21,-5,20])
plt.title("KNeighborsRegressor (k = %i)" % (n_neighbors))
plt.xlabel('Target Velocity (degree/s)')
plt.ylabel('Eye Velocity (degree/s)')
plt.show()



###############################################################################
# knearest regression x-validation
n_neighbors = 10
kf = KFold(80, n_folds=80)
Error = np.empty([80,1])
i = 0
for train, test in kf:
    
    X_train, X_test, y_train, y_test = Xd1[train], Xd1[test], Yd1[train], Yd1[test]
    
    knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
    y1_ = knn.fit(X_train, y_train).predict(X_test)
    Error[i] = (y1_-y_test)**2
    i += 1
    
Error1 = np.sqrt(np.mean(Error))

i = 0
for train, test in kf:
    
    X_train, X_test, y_train, y_test = Xd6[train], Xd6[test], Yd6[train], Yd6[test]
    
    knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
    y1_ = knn.fit(X_train, y_train).predict(X_test)
    Error[i] = (y1_-y_test)**2
    i += 1
    
Error6 = np.sqrt(np.mean(Error))

i = 0
for train, test in kf:
    
    X_train, X_test, y_train, y_test = Xd10[train], Xd10[test], Yd10[train], Yd10[test]
    
    knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
    y1_ = knn.fit(X_train, y_train).predict(X_test)
    Error[i] = (y1_-y_test)**2
    i += 1
    
Error10 = np.sqrt(np.mean(Error))

i = 0
for train, test in kf:
    
    X_train, X_test, y_train, y_test = Xd20[train], Xd20[test], Yd20[train], Yd20[test]
    
    knn = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
    y1_ = knn.fit(X_train, y_train).predict(X_test)
    Error[i] = (y1_-y_test)**2
    i += 1
    
Error20 = np.sqrt(np.mean(Error))
fig = plt.figure()
plt.plot([1,6,10,20], [Error1,Error6,Error10,Error20], '-k', [1,6,10,20], [Error1,Error6,Error10,Error20], 'ob')
plt.xlabel('Dots Diameter');plt.ylabel('X-validation MSE');plt.title('KNeighborsRegressor')
plt.show()

###############################################################################
# Bayesian Regression
clf = linear_model.BayesianRidge(compute_score=True,normalize=True)
#y1_ = clf.fit(Xd1, Yd1).predict(T)
#c1 = clf.coef_
#a1 = clf.alpha_
#l1 = clf.lambda_

y6_ = clf.fit(Xd6, Yd6).predict(T)
c6 = clf.coef_
a6 = clf.alpha_
l6 = clf.lambda_

y10_ = clf.fit(Xd10, Yd10).predict(T)
c10 = clf.coef_
a10 = clf.alpha_
l10 = clf.lambda_

y20_ = clf.fit(Xd20, Yd20).predict(T)
c20 = clf.coef_
a20 = clf.alpha_
l20 = clf.lambda_

#fig1 = plt.figure()
#plt.scatter(Xd1, Yd1, c='k', label='data')
#line_1, = plt.plot(T, y1_, 'k', label='1 degree')
#plt.text(10, 18, r'$p(V_e|V_t,\omega,\alpha) = \mathcal{N}(V_e|V_t\omega,\alpha)$')
#plt.text(10, 17, r'$p(\omega|\lambda) = \mathcal{N}(\omega|0,\lambda^{-1})$')
#plt.text(10, 16, r'$\omega = %.3f$' % (c1))
#plt.text(10, 15, r'$\alpha = %.3f$' % (a1))
#plt.text(10, 14, r'$\lambda = %.3f$' % (l1))
#plt.axis([9,21,-5,20])
#plt.title("Bayesian Regression - d = 1 degree")
#plt.xlabel('Target Velocity (degree/s)')
#plt.ylabel('Eye Velocity (degree/s)')
#
#plt.show()

fig2 = plt.figure()
plt.scatter(Xd6, Yd6, c='b', label='data')
line_1, = plt.plot(T, y6_, 'b', label='1 degree')
plt.text(2, 18, r'$p(V_e|V_t,\omega,\alpha) = \mathcal{N}(V_e|V_t\omega,\alpha)$')
plt.text(2, 17, r'$p(\omega|\lambda) = \mathcal{N}(\omega|0,\lambda^{-1})$')
plt.text(2, 16, r'$\omega = %.3f$' % (c6))
plt.text(2, 15, r'$\alpha = %.3f$' % (a6))
plt.text(2, 14, r'$\lambda = %.3f$' % (l6))
plt.axis([0,25,-5,20])
plt.title("Bayesian Regression - 6 degrees")
plt.xlabel('Target Velocity (degree/s)')
plt.ylabel('Eye Velocity (degree/s)')
plt.show()

fig3 = plt.figure()
plt.scatter(Xd10, Yd10, c='g', label='data')
line_1, = plt.plot(T, y10_, 'g', label='1 degree')
plt.text(2, 18, r'$p(V_e|V_t,\omega,\alpha) = \mathcal{N}(V_e|V_t\omega,\alpha)$')
plt.text(2, 17, r'$p(\omega|\lambda) = \mathcal{N}(\omega|0,\lambda^{-1})$')
plt.text(2, 16, r'$\omega = %.3f$' % (c10))
plt.text(2, 15, r'$\alpha = %.3f$' % (a10))
plt.text(2, 14, r'$\lambda = %.3f$' % (l10))
plt.axis([0,25,-5,20])
plt.title("Bayesian Regression - d = 10 degrees")
plt.xlabel('Target Velocity (degree/s)')
plt.ylabel('Eye Velocity (degree/s)')
plt.show()

fig4 = plt.figure()
plt.scatter(Xd20, Yd20, c='r', label='data')
line_1, = plt.plot(T, y20_, 'r', label='1 degree')
plt.text(2, 18, r'$p(V_e|V_t,\omega,\alpha) = \mathcal{N}(V_e|V_t\omega,\alpha)$')
plt.text(2, 17, r'$p(\omega|\lambda) = \mathcal{N}(\omega|0,\lambda^{-1})$')
plt.text(2, 16, r'$\omega = %.3f$' % (c20))
plt.text(2, 15, r'$\alpha = %.3f$' % (a20))
plt.text(2, 14, r'$\lambda = %.3f$' % (l20))
plt.axis([0,25,-5,20])
plt.title("Bayesian Regression - d = 20 degrees")
plt.xlabel('Target Velocity (degree/s)')
plt.ylabel('Eye Velocity (degree/s)')
plt.show()

###############################################################################
# Bayesian Regression x-validation
NumSamples = 206;
clf = linear_model.BayesianRidge(compute_score=True)
kf = KFold(NumSamples, n_folds=NumSamples)
c = np.empty([NumSamples,1])
a = np.empty([NumSamples,1])
l = np.empty([NumSamples,1])
Error = np.empty([NumSamples,1])
i = 0
#for train, test in kf:
#    
#    X_train, X_test, y_train, y_test = Xd1[train], Xd1[test], Yd1[train], Yd1[test]
#    clf.fit(X_train, y_train)
#    y1_ = clf.predict(X_test)
#    Error[i] = (y1_-y_test)**2
#    c[i] = clf.coef_
#    a[i] = clf.alpha_
#    l[i] = clf.lambda_
#    i += 1
#    
#Error1 = np.sqrt(np.mean(Error))
#c1 = np.sqrt(np.mean(c))
#a1 = np.sqrt(np.mean(a))
#l1 = np.sqrt(np.mean(l))


i = 0
for train, test in kf:
    
    X_train, X_test, y_train, y_test = Xd6[train], Xd6[test], Yd6[train], Yd6[test]
    clf.fit(X_train, y_train)
    y1_ = clf.predict(X_test)
    Error[i] = (y1_-y_test)**2
    c[i] = clf.coef_
    a[i] = clf.alpha_
    l[i] = clf.lambda_
    i += 1
    
Error6 = np.sqrt(np.mean(Error))
c6 = (np.nanmean(c))
a6 = (np.median(a))
l6 = (np.median(l))
Error6s = np.sqrt(np.std(Error))
c6s = (np.std(c))
a6s = (np.std(a))
l6s = (np.std(l))

kf = KFold(NumSamples, n_folds=NumSamples)
c = np.empty([NumSamples,1])
a = np.empty([NumSamples,1])
l = np.empty([NumSamples,1])
Error = np.empty([NumSamples,1])
i = 0
for train, test in kf:
    
    X_train, X_test, y_train, y_test = Xd10[train], Xd10[test], Yd10[train], Yd10[test]
    clf.fit(X_train, y_train)
    y1_ = clf.predict(X_test)
    Error[i] = (y1_-y_test)**2
    c[i] = clf.coef_
    a[i] = clf.alpha_
    l[i] = clf.lambda_
    i += 1
    
Error10 = np.sqrt(np.mean(Error))
#Error10 = Error
c10 = (np.nanmean(c))
a10 = (np.median(a))
l10 = (np.median(l))
Error10s = np.sqrt(np.std(Error))
c10s = (np.std(c))
a10s = (np.std(a))
l10s = (np.std(l))

kf = KFold(NumSamples, n_folds=NumSamples)
c = np.empty([NumSamples,1])
a = np.empty([NumSamples,1])
l = np.empty([NumSamples,1])
Error = np.empty([NumSamples,1])
i = 0
for train, test in kf:
    
    X_train, X_test, y_train, y_test = Xd20[train], Xd20[test], Yd20[train], Yd20[test]
    clf.fit(X_train, y_train)
    y1_ = clf.predict(X_test)
    Error[i] = (y1_-y_test)**2
    c[i] = clf.coef_
    a[i] = clf.alpha_
    l[i] = clf.lambda_
    i += 1
    
Error20 = np.sqrt(np.mean(Error))
c20 = (np.nanmean(c))
a20 = (np.median(a))
l20 = (np.median(l))
Error20s = np.sqrt(np.std(Error))
c20s = (np.std(c))
a20s = (np.std(a))
l20s = (np.std(l))

fig = plt.figure()
plt.plot([6,10,20], [Error6,Error10,Error20], '-k', [6,10,20], [Error6,Error10,Error20], 'ob')
plt.xlabel('Dots Diameter');plt.ylabel('X-validation MSE');plt.title('Bayesian Regression')
plt.show()

fig2 = plt.figure()
plt.plot([6,10,20], [c6,c10,c20], '-k', [6,10,20], [c6,c10,c20], 'ob')
plt.xlabel('Dots Diameter');plt.ylabel('coefficient');plt.title('Bayesian Regression')
plt.grid(True)
plt.show()

fig3 = plt.figure()
plt.plot([6,10,20], [a6,a10,a20], '-k', [6,10,20], [a6,a10,a20], 'ob')
plt.xlabel('Dots Diameter');plt.ylabel('alpha');plt.title('Bayesian Regression')
plt.grid(True)
plt.show()

fig4 = plt.figure()
plt.plot([6,10,20], [l6,l10,l20], '-k', [6,10,20], [l6,l10,l20], 'ob')
plt.xlabel('Dots Diameter');plt.ylabel('lambda');plt.title('Bayesian Regression')
plt.grid(True)
plt.show()

###############################################################################
# Multi-subject figures
#Error6_sb = Error6
#Error10_sb = Error10
#Error20_sb = Error20
#c6_sb = c6
#c10_sb = c10
#c20_sb = c20
#a6_sb = a6
#a10_sb = a10
#a20_sb = a20
#
#Error6_az = Error6
#Error10_az = Error10
#Error20_az = Error20
#c6_az = c6
#c10_az = c10
#c20_az = c20
#a6_az = a6
#a10_az = a10
#a20_az = a20
#
#Error6_gc = Error6
#Error10_gc = Error10
#Error20_gc = Error20
#c6_gc = c6
#c10_gc = c10
#c20_gc = c20
#a6_gc = a6
#a10_gc = a10
#a20_gc = a20
#
#fig = plt.figure()
#plt.plot([2,6,20], [Error6_sb,Error10_sb,Error20_sb], '-k', [6,10,20], [Error6_sb,Error10_sb,Error20_sb], 'ob',
#        [2,6,20], [Error6_az,Error10_az,Error20_az], '-k', [6,10,20], [Error6_az,Error10_az,Error20_az], 'ob',
#        [2,6,20], [Error6_gc,Error10_gc,Error20_gc], '-k', [6,10,20], [Error6_gc,Error10_gc,Error20_gc], 'ob'
#        )
#plt.xlabel('Dots Diameter');plt.ylabel('X-validation MSE');plt.title('Bayesian Regression')
#plt.show()
#
#fig2 = plt.figure()
#plt.plot([6,10,20], [c6_sb,c10_sb,c20_sb], '-k', [6,10,20], [c6_sb,c10_sb,c20_sb], 'ob',
#        [6,10,20], [c6_az,c10_az,c20_az], '-k', [6,10,20], [c6_az,c10_az,c20_az], 'ob',
#        [6,10,20], [c6_gc,c10_gc,c20_gc], '-k', [6,10,20], [c6_gc,c10_gc,c20_gc], 'ob'
#        )
#plt.xlabel('Dots Diameter');plt.ylabel('coefficient');plt.title('Bayesian Regression')
#plt.grid(True)
#plt.show()
#
#fig3 = plt.figure()
#plt.plot([6,10,20], [a6_sb,a10_sb,a20_sb], '-k', [6,10,20], [a6_sb,a10_sb,a20_sb], 'ob',
#        [6,10,20], [a6_az,a10_az,a20_az], '-k', [6,10,20], [a6_az,a10_az,a20_az], 'ob',
#        [6,10,20], [a6_gc,a10_gc,a20_gc], '-k', [6,10,20], [a6_gc,a10_gc,a20_gc], 'ob'
#        )
#plt.xlabel('Dots Diameter');plt.ylabel('alpha');plt.title('Bayesian Regression')
#plt.grid(True)
#plt.show()
#
#
#
#
