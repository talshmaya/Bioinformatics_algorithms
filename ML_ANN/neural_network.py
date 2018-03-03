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
from sklearn.feature_selection import SelectFromModel
import sys
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.feature_selection import SelectPercentile
import os
from sklearn.ensemble import ExtraTreesClassifier
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import sklearn
print sklearn.__version__
from sklearn.model_selection import LeaveOneOut

def run_NN():
	# fix random seed for reproducibility
	seed = 7
	numpy.random.seed(seed)
	# load dataset
	file_name=sys.argv[1]
	dataset= numpy.loadtxt(file_name,delimiter=",")
	# split into input (X) and output (Y) variables

	#X = dataset[0:-160,0:-2]
	#Y = dataset[0:-160,-1]
	X = dataset[:,0:-2]
	Y = dataset[:,-1]

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

	for j in range(0,1):
		co=0
		X = dataset[101:-1,0:-2]
		Y = dataset[101:-1,-1]
	#	test = SelectKBest(score_func=SelectPercentile) #, k=int(NUM_FEAT)
		test = SelectKBest(score_func=chi2, k=int(NUM_FEAT))
		fit = test.fit(X, Y)


		get_features(test)

		# summarize scores
		numpy.set_printoptions(precision=3)
		variables=[]
		old_j=0
		features = fit.transform(X)
		# summarize selected features
		#print(numpy.transpose(features).shape)
		#X_new=features

		X = features
		DIM=len(X[1])
		loo = LeaveOneOut()
		tp_a, tn_a , fp_a, fn_a=[],[],[],[]
		#for train, test in loo.split(X):

		# create model
		model = Sequential()
		model.add(Dense(int(NUM_NEU), activation="relu", kernel_initializer="uniform", input_dim=DIM))
		#model.add(Dense(12, init='uniform', activation='relu'))
		model.add(Dense(1, init='uniform', activation='sigmoid'))
		# Compile model
		model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
		# Fit the model

		# split into input (X) and output (Y) variables
		X_test = dataset[0:100,0:-2]
		Y_test = dataset[0:100,-1]

		model.fit(X, Y, nb_epoch=int(NUM_EP), batch_size=10,  verbose=0)
		# calculate predictions


		#X_test = dataset[0:99,0:-2]

		b=fit.transform(X_test)
		#Y_test = dataset[test,-1]
		predictions = model.predict(b)
		# round predictions
		rounded = [round(x[0]) for x in predictions]
		print ("prediction is: ")
		#for i,x in enumerate(predictions):
		#	print x, Y[i]
		p=[int(x) for x in rounded]
		p=''.join(map(str,p))
		print p
		#rounded_test_label = [round(x[0]) for x in Y_test]
		print ("observation is: ")

		o=[int(x) for x in Y_test]
		o=''.join(map(str,o))
		print o
		import difflib
		fp,fn,tp,tn=0,0,0,0
		print ("----------------")
		co+=1
		for O,P in zip(o,p):
			if O=="1" and P=="1":
				tp+=1
				tp_a.append(True)
			if  O=='0' and P=='0':
				tn+=1
				tn_a.append(True)
			if O=='1' and P=='0':
				fp+=1
				fp_a.append(True)
			if O=='0' and P=='1': # prediction is yes but ovservation is no
				fn+=1
				fn_a.append(True)

		print ("false positive and false negatives: ")
		print  fp, fn
		print("true positives and true negatives: ")
		print tp,tn

		scores = model.evaluate(X,Y)
		print ("----------------")
		print "1 layer with " + str(NUM_NEU) + " nuerons. " + str(NUM_EP) + " Epochs. " + str(NUM_FEAT) + " features were selected. "
		print("Training Accuracy:")
		#print("\n%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))
		test_score=float(tp+tn)/len(Y_test)*100
		test_all.append(test_score)
		train_score=scores[1]*100
		train_all.append(train_score)
		print(scores[1]*100)
		print "Testing Accuracy: "
		print (float(tp+tn)/len(Y_test)*100)
	"""
	print ("+++++++++++++++++++++++++++++++++")
	print (len(tp_a))+(len(tn_a))
	print (len(fp_a))+(len(fn_a))
	print "Total Accuracy: "
	print float( ( (len(tp_a))+(len(tn_a)) ))/ co
	#print "Max training score for " + str(NUM_FEAT) + " features: " + str(max(train_all))
	#print "Max testing score for " + str(NUM_FEAT) + " features: " + str(max(test_all))
	train_all,test_all=[],[]
	train_score,test_score=0,0
	"""
def get_features(t):
	mutation_file = open('/home/tshmaya/Workspace/ML_drug_resistance/rif_unexplained_muts_new')
	for i in mutation_file:
		mut=i.split('\t')
	supporting_vars=t.get_support()
	print supporting_vars
	print len(supporting_vars)
	print mut[22563]
	print mut[2]
	#for inx,i in enumerate(supporting_vars):
		#if i:
			#print inx
			#print mut[inx]


if __name__ == "__main__":

	run_NN()