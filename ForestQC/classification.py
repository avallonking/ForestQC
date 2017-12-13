# train a random forest classifier to identify good or bad variants from grey variants

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support
import pandas as pd
import os

# just a second choice
def random_forest_classifierA(labelled_data, grey_variants, user_features):
    # input: a dataset that has balanced sample size of good and bad variants
    # output: predicted good or bad variants from grey variants
    try:
        _features = ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','ABHet', 'ABHom', 'GC'] + \
                    user_features
    except TypeError:
        _features = ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','ABHet', 'ABHom', 'GC']
    x, y = labelled_data.loc[:, _features].values, labelled_data.loc[:, 'Good'].values
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=1)
    rf = RandomForestClassifier(random_state=1, n_jobs=8, n_estimators=50)

    print('\nTraining model...')
    rf.fit(x_train, y_train)

    weights = rf.feature_importances_
    print('Feature','Weight')
    for feature, weight in zip(_features, weights):
      print(feature, weight)

    print('\nTesting model...')
    precison_recall = precision_recall_fscore_support(y_true=y_test, y_pred=rf.predict(x_test))
    print('\n\t\tBad\tGood')
    print('Precision: ' + str(precison_recall[0]))
    print('Recall: ' + str(precison_recall[1]))
    print('F1-score: ' + str(precison_recall[2]))
    print('Sample size: ' + str(precison_recall[3]))
    
    x_pred = grey_variants.loc[:, _features].values
    print('\nPredicting variants...')
    y_pred = rf.predict(x_pred)
    print('Done.\n')

    return y_pred

def random_forest_classifierB(labelled_data, grey_variants, user_features):
    # input: a dataset that has balanced sample size of good and bad variants
    # output: predicted good or bad variants from grey variants

    try:
        _features = ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','GC'] + user_features
    except TypeError:
        _features = ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','GC']

    x, y = labelled_data.loc[:, _features].values, labelled_data.loc[:, 'Good'].values
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=1)
    rf = RandomForestClassifier(random_state=1, n_jobs=8, n_estimators=50)
    print('\nTraining model...')
    rf.fit(x_train, y_train)

    weights = rf.feature_importances_
    print('Feature','Weight')
    for feature, weight in zip(_features, weights):
      print(feature, weight)

    print('\nTesting model...')
    precison_recall = precision_recall_fscore_support(y_true=y_test, y_pred=rf.predict(x_test))
    print('\n\t\tBad\tGood')
    print('Precision: ' + str(precison_recall[0]))
    print('Recall: ' + str(precison_recall[1]))
    print('F1-score: ' + str(precison_recall[2]))
    print('Sample size: ' + str(precison_recall[3]))

    x_pred = grey_variants.loc[:, _features].values
    print('\nPredicting variants...')
    y_pred = rf.predict(x_pred)
    print('Done.')

    return y_pred

def classification(good, bad, grey, model, user_features):
  rf_model = {'A': random_forest_classifierA, 'B': random_forest_classifierB}
  if good.shape[0] > bad.shape[0]:
    prediction = rf_model[model](pd.concat([good.sample(n=bad.shape[0],random_state=9), bad]), grey, user_features)
  elif good.shape[0] == bad.shape[0]:
    prediction = rf_model[model](pd.concat([good, bad]), grey, user_features)
  else:
    prediction = rf_model[model](pd.concat([bad.sample(n=good.shape[0],random_state=9), good]), grey, user_features)

  return prediction


def execute_classification(good_variants, bad_variants, grey_variants, model, output_handle, user_features):
  dir = os.path.dirname(os.path.realpath(good_variants))
  print('Loading data...')
  good = pd.read_table(good_variants)
  bad = pd.read_table(bad_variants)
  grey = pd.read_table(grey_variants)

  pred = classification(good, bad, grey, model, user_features)
  grey['Good'] = pred
  
  predicted_good = grey[grey['Good'] == 1]
  predicted_bad = grey[grey['Good'] == 0]
  
  print('Number of predicted good variants: ' + str(predicted_good.shape[0]))
  print('Number of predicted bad variants: ' + str(predicted_bad.shape[0]))

  print('\nWriting data...')
  predicted_good.to_csv(dir + '/' + 'predicted_good.' + output_handle, header=None, index = False, sep = '\t', na_rep='NA')
  predicted_bad.to_csv(dir + '/' + 'predicted_bad.' + output_handle, header=None, index = False, sep = '\t', na_rep='NA')
  print('Done.')
