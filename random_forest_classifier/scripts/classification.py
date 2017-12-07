# train a random forest classifier to identify good or bad variants from grey variants

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support
import pandas as pd
import os

def random_forest_classifierB(labelled_data, grey_variants):
    # input: a dataset that has balanced sample size of good and bad variants
    # output: predicted good or bad variants from grey variants
    x, y = labelled_data.loc[:, ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','ABHet', 'ABHom', 'GC']].values, labelled_data.loc[:, 'Good'].values
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=1)
    rf = RandomForestClassifier(random_state=1, n_jobs=8, n_estimators=50)

    print('\nTraining model...')
    rf.fit(x_train, y_train)
    print('Done.\n')

    weights = rf.feature_importances_
    features = ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','ABHet', 'ABHom', 'GC']
    print('Feature','Weight')
    for feature, weight in zip(features, weights):
      print(feature, weight)

    print('\nTesting model...')
    precison_recall = precision_recall_fscore_support(y_true=y_test, y_pred=rf.predict(x_test))
    print('Done')
    print('\n\t\tBad\tGood')
    print('Precision: ' + str(precison_recall[0]))
    print('Recall: ' + str(precison_recall[1]))
    print('F1-score: ' + str(precison_recall[2]))
    print('Sample size: ' + str(precison_recall[3]))
    
    x_pred = grey_variants.loc[:, ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','ABHet', 'ABHom', 'GC']].values
    print('\nPredicting variants...')
    y_pred = rf.predict(x_pred)
    print('Done.\n')

    return y_pred

# deprecated
def random_forest_classifierA(labelled_data, grey_variants):
    # input: a dataset that has balanced sample size of good and bad variants
    # output: predicted good or bad variants from grey variants
    x, y = labelled_data.loc[:, ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','GC']].values, labelled_data.loc[:, 'Good'].values
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=1)
    # pipeline = Pipeline([('scl', StandardScaler()), ('rf', RandomForestClassifier(random_state=1,n_jobs=8,n_estimators=10))])
    rf = RandomForestClassifier(random_state=1, n_jobs=8, n_estimators=50)
    print('\nTraining model...')
    rf.fit(x_train, y_train)
    print('Done.\n')

    weights = rf.feature_importances_
    features = ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','GC']
    print('Feature','Weight')
    for feature, weight in zip(features, weights):
      print(feature, weight)

    print('\nTesting model...')
    precison_recall = precision_recall_fscore_support(y_true=y_test, y_pred=rf.predict(x_test))
    print('Done')
    print('\n\t\tBad\tGood')
    print('Precision: ' + str(precison_recall[0]))
    print('Recall: ' + str(precison_recall[1]))
    print('F1-score: ' + str(precison_recall[2]))
    print('Sample size: ' + str(precison_recall[3]))

    x_pred = grey_variants.loc[:, ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','GC']].values
    print('\nPredicting variants...')
    y_pred = rf.predict(x_pred)
    print('Done.')

    return y_pred

def classification(good, bad, grey, model = 'A'):
  rf_model = {'A': random_forest_classifierA, 'B': random_forest_classifierB}
  if good.shape[0] > bad.shape[0]:
    prediction = rf_model[model](pd.concat([good.sample(n=bad.shape[0],random_state=9), bad]), grey)
  elif good.shape[0] == bad.shape[0]:
    prediction = rf_model[model](pd.concat([good, bad]), grey)
  else:
    prediction = rf_model[model](pd.concat([bad.sample(n=good.shape[0],random_state=9), good]), grey)

  return prediction


def executeClassification(good_variants, bad_variants, grey_variants, model = 'A', output_handle = 'variants.tsv'):
  dir = os.path.dirname(good_variants)
  print('Loading data...')
  good = pd.read_table(good_variants, header=None)
  bad = pd.read_table(bad_variants, header=None)
  grey = pd.read_table(grey_variants, header=None)
  columns1 =  ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom', 'GC', 'Good']
  columns2 =  ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom', 'GC']
  good.columns = columns1
  bad.columns = columns1
  grey.columns = columns2
  print('Done.')

  pred = classification(good, bad, grey, model)
  grey['Good'] = pred
  
  predicted_good = grey[grey['Good'] == 1]
  predicted_bad = grey[grey['Good'] == 0]
  
  print('Number of predicted good variants: ' + str(predicted_good.shape[0]))
  print('Number of predicted bad variants: ' + str(predicted_bad.shape[0]))

  print('\nWriting data...')
  predicted_good.to_csv(dir + '/' + 'predicted_good.' + output_handle, index = False, sep = '\t', na_rep='NA')
  predicted_bad.to_csv(dir + '/' + 'predicted_bad.' + output_handle, index = False, sep = '\t', na_rep='NA')
  print('Done.')
