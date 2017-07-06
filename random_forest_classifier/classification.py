# train a random forest classifier to identify good or bad variants from grey variants

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_fscore_support
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import pandas as pd
import sys

def random_forest_classifier(labelled_data, grey_variants):
    # input: a dataset that has balanced sample size of good and bad variants
    # output: predicted good or bad variants from grey variants
    x, y = labelled_data.loc[:, ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','ABHet']].values, labelled_data.loc[:, 'Good'].values
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=1)
    pipeline = Pipeline([('scl', StandardScaler()), ('rf', RandomForestClassifier(random_state=1,n_jobs=8,n_estimators=10))])
    pipeline.fit(x_train, y_train)
    precison_recall = precision_recall_fscore_support(y_true=y_test, y_pred=pipeline.predict(x_test))
    print('\tBad\tGood')
    print('Precision: ' + str(precison_recall[0]))
    print('Recall: ' + str(precison_recall[1]))
    print('F1-score: ' + str(precison_recall[2]))
    print('Sample size: ' + str(precison_recall[3]))

    x_pred = grey_variants.loc[:, ['Mean_DP','Mean_GQ','SD_DP','SD_GQ','Outlier_DP','Outlier_GQ','ABHet']].values
    y_pred = pipeline.predict(x_pred)

    return y_pred

def main():
  good_variants = sys.argv[1]
  bad_variants = sys.argv[2]
  grey_variants = sys.argv[3]
  output_handle = sys.argv[4] if len(sys.argv) >= 5 else 'variants.tsv'

  good = pd.read_table(good_variants, header=None)
  bad = pd.read_table(bad_variants, header=None)
  grey = pd.read_table(grey_variants, header=None)

  columns1 =  ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom', 'Good']
  columns2 =  ['RSID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'Mean_DP', 'Mean_GQ', 'SD_DP', 'SD_GQ', 'Outlier_DP', 'Outlier_GQ', 'Discordant_Geno', 'Mendel_Error', 'Missing_Rate', 'HWE', 'ABHet', 'ABHom']
  good.columns = columns1
  bad.columns = columns1
  grey.columns = columns2

  if good.shape[0] > bad.shape[0]:
    prediction = random_forest_classifier(pd.concat([good.sample(n=bad.shape[0]), bad]), grey)
  elif good.shape[0] == bad.shape:
    prediction = random_forest_classifier(pd.concat([good, bad]), grey)
  else:
    prediction = random_forest_classifier(pd.concat([bad.sample(n=good.shape[0]), good]), grey)

  grey['Good'] = prediction
  predicted_good = grey[grey['Good'] == 1]
  predicted_bad = grey[grey['Good'] == 0]
  print('Predicted good variants: ' + str(predicted_good.shape[0]))
  print('Predicted bad variants: ' + str(predicted_bad.shape[0]))
  predicted_good.to_csv('predicted_good_' + output_handle, index = False, sep = '\t', na_rep='NA')
  predicted_bad.to_csv('predicted_bad_' + output_handle, index = False, sep = '\t', na_rep='NA')

if __name__ == '__main__':
  main()
