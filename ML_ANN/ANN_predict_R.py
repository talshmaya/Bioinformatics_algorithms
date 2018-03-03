# Tal Shmaya
# to run:
# $ python NN.py input_file num_features num_neurons num_epochs
# example:
# tshmaya@valafar09:~/Workspace/ML_drug_resistance$ python NN.py "shuf_inh_equal" 10 1 50

from keras.models import Sequential
from keras.layers import Dense
import numpy
import difflib
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
#from sklearn.datasets import load_iris
#from sklearn.feature_selection import SelectFromModel
import sys
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
import os
from sklearn.ensemble import ExtraTreesClassifier
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from sklearn.model_selection import LeaveOneOut
# fix random seed for reproducibility
seed = 7
numpy.random.seed(seed)
# load pima indians dataset
file_name=sys.argv[1]
#dataset = numpy.loadtxt("pima-indians-diabetes-head.csv", delimiter=",")
dataset= numpy.loadtxt(file_name,delimiter=",")
# split into input (X) and output (Y) variables

#X = dataset[0:-160,0:-2]
#Y = dataset[0:-160,-1]
X = dataset[100:-1,0:-2]
Y = dataset[100:-1,-1]

NUM_FEAT=sys.argv[2]
NUM_NEU=sys.argv[3]
NUM_EP=sys.argv[4]


"""
model = ExtraTreesClassifier()
model.fit(X, Y)
print([i for i in model.feature_importances_])
print sorted(model.feature_importances_)
print max(i for i in model.feature_importances_)


print X.shape
#lsvc = LinearSVC(C=0.01, penalty="l1", dual=False).fit(X, Y)
#model = SelectFromModel(lsvc, prefit=True)
X_new = SelectKBest(chi2, k=2).fit(X, Y)
#X_new = model.transform(X)
print ([i for i in X_new])
#DIM=len(X_new)-2
"""
train_all,test_all=[],[]
train_score,test_score=0,0

# feature extraction
#feat=[1,2,3,4,5,6,7,8,9,10,15,20,50,100,150,200,300,355]
feat=[1000,1300,1500,2000,2500,3000]
for i in feat:
	for j in range(0,4):
		X = dataset[100:-1,0:-2]
		Y = dataset[100:-1,-1]
		test = SelectKBest(score_func=chi2, k=i)
		fit = test.fit(X, Y)
		# summarize scores
		numpy.set_printoptions(precision=3)
		variables=[]
		old_j=0

		"""
		for i,j in enumerate(fit.scores_):
				variables.append((i,float(j)))
				if j > 60:
					print i+1
					max_j=j
				j=old_j

		#print(sorted(variables,key=lambda student: float(student[1])))
		"""
		features = fit.transform(X)
		# summarize selected features
		#print(numpy.transpose(features).shape)
		#X_new=features

		X = features

		DIM=len(X[1])

		"""
		model = LogisticRegression()
		rfe = RFE(model, 3)
		fit = rfe.fit(X, Y)
		print("Num Features: %d") % fit.n_features_
		print("Selected Features: %s") % fit.support_
		print("Feature Ranking: %s") % fit.ranking_
		"""

		# create model
		model = Sequential()
		#model.add(Dense(int(NUM_NEU), input_dim=DIM, init='uniform', activation='relu'))
		model.add(Dense(int(i), activation="relu", kernel_initializer="uniform", input_dim=DIM))
		#model.add(Dense(12, init='uniform', activation='relu'))
		#model.add(Dense(12, init='uniform', activation='relu'))
		#model.add(Dense(1, init='uniform', activation='sigmoid'))
		model.add(Dense(1, activation="sigmoid", kernel_initializer="uniform"))
		# Compile model
		model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
		# Fit the model

		# split into input (X) and output (Y) variables
		#X_test = dataset[-161:-1,0:-2]
		#Y_test = dataset[-161:-1,-1]



		model.fit(X, Y, nb_epoch=int(NUM_EP), batch_size=10,  verbose=0)
		# calculate predictions


		X_test = dataset[0:99,0:-2]

		b=fit.transform(X_test)
		Y_test = dataset[0:99,-1]
		predictions = model.predict(b)
		# round predictions
		rounded = [round(x[0]) for x in predictions]
		#print ("prediction is: ")
		p=[int(x) for x in rounded]
		p=''.join(map(str,p))
		#print p
		#rounded_test_label = [round(x[0]) for x in Y_test]
		#print ("observation is: ")

		o=[int(x) for x in Y_test]
		o=''.join(map(str,o))
		#print o
		#import difflib
		fp,fn,tp,tn=0,0,0,0
		for O,P in zip(o,p):
			if O=="1" and P=="1":
				tp+=1
			if  O=='0' and P=='0':
				tn+=1
			if O=='1' and P=='0':
				fp+=1
			if O=='0' and P=='1': # prediction is yes but ovservation is no
				fn+=1
		"""
		print ("false positive and false negatives: ")
		print  fp, fn
		print("true positives and true negatives: ")
		print tp,tn
		"""
		scores = model.evaluate(X ,Y)
		print ("----------------")
		print "1 layer with " + str(NUM_NEU) + " nuerons. " + str(NUM_EP) + " Epochs. " + str(i) + " features were selected. "
		print("Training Accuracy:")
		#print("\n%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))
		test_score=float(tp+tn)/len(Y_test)*100
		test_all.append(test_score)
		train_score=scores[1]*100
		train_all.append(train_score)
		print(scores[1]*100)
		print "Testing Accuracy: "
		print (float(tp+tn)/len(Y_test)*100)
	print "Max training score for " + str(i) + " features: " + str(max(train_all))
	print "Max testing score for " + str(i) + " features: " + str(max(test_all))
	train_all,test_all=[],[]
	train_score,test_score=0,0