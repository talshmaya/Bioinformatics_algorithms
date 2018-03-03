import numpy as np
from sklearn.ensemble import RandomForestClassifier

#initialize data set and classification
X,Y=[],[]

#get current data
with open ('current_data.csv','r') as f:
	lines=f.readlines()
for line in lines:	
		l=line.split(',')
		X.append(l)

#get classification
with open ('current_classifications.csv','r') as f:
	lines=f.readlines()
for line in lines:
	Y=line.split(',')

#construct the forest with 200 trees
clf = RandomForestClassifier(n_estimators=200)
clf = clf.fit(X, Y)

#get importance of features
importances = clf.feature_importances_
indices = np.argsort(importances)[::-1]
print [i for i in indices]
#print  clf.score(X, Y)

# classify a new dataset using the model we trained above
X_test=[]
with open ('new_data.csv','r') as f:
	lines=f.readlines()
for line in lines:	
		l=line.split(',')
		X_test.append(l)

print clf.predict(X_test)