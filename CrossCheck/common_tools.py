import pandas as pd
import numpy as np
import json
from sklearn.metrics import roc_auc_score, f1_score, roc_curve, average_precision_score
from pprint import pprint 
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from matplotlib.patches import Patch, Circle # Correcting the patch error
from matplotlib.lines import Line2D
import uproot3 as uproot

def correlation_to_covariance(correlation_matrix, std_deviations):
    # Calculate the covariance matrix
    cov_matrix = correlation_matrix * np.outer(std_deviations, std_deviations)
    
    return cov_matrix

def create_dict_from_df_destructive(dataframe, destructive=True):
    dict_df = dict()
    dataframe.reset_index(inplace=True)
    if 'Slice' in dataframe:
        for trg in set(dataframe.Slice):
            print(f'Slice_HLT_{trg}')
            dataframe[f'Slice_HLT_{trg}'] = 0
            dataframe[f'Slice_HLT_{trg}'][dataframe.Slice==trg] = 1
        dataframe.drop('Slice', axis=1, inplace=True)
        
    #pdb.set_trace()
    for var in dataframe:
        if 'HLT' in var: dict_df[var] = np.array(dataframe[var], dtype=np.int8)
        if var in ['run', 'luminosityBlock', 'event', 'subentry']:
            print(var)
            dict_df[var] = np.array(dataframe[var], dtype=np.int)
        else: dict_df[var] = np.array(dataframe[var])
        if destructive: dataframe.drop(var, axis=1, inplace=True)
    return dict_df


def save_dict_to_root(dictionary, output_file, tree):
    dict_types = dict()
    for var in dictionary: 
        if 'object' in str(type(dictionary[var].dtype)):
            print(var)
            dictionary[var] = list(dictionary[var])
            dict_types[var] = str
        else: dict_types[var] = dictionary[var].dtype
    uproot_file = uproot.recreate(output_file)
    uproot_file[tree] = uproot.newtree(dict_types)
    uproot_file[tree].extend(dictionary)
    uproot_file.close()


def save_df_to_root(dataframe, output_file, tree):
    dictionary = create_dict_from_df_destructive(dataframe)
    save_dict_to_root(dictionary, output_file, tree)


def save_data(dataframe, path, **kwargs):
    if path.endswith('csv'):
        dataframe.to_csv(path)
    elif path.endswith('pkl'):
        dataframe.to_pickle(path)
    elif path.endswith('root'):
        save_df_to_root(dataframe, path, **kwargs)





def params_to_string(model):
  model.get_params()
  string = ''
  lens = [len(k) for k in model.get_params()]
  max_len = max(lens)
  for k,v in model.get_params().items():
    #string += f'{k+" "*(max_len-len(k))} =  {v}\n'
    string += f'{k} = {v}\n'

  return string


def plot_classifier_distributions(model, test, train, cols, print_params=False, params=None):

    test_background = model.predict_proba(test.query('label==0')[cols])[:,1]
    test_signal     = model.predict_proba(test.query('label==1')[cols])[:,1]
    train_background= model.predict_proba(train.query('label==0')[cols])[:,1]
    train_signal    = model.predict_proba(train.query('label==1')[cols])[:,1]

    test_pred = model.predict_proba(test[cols])[:,1]
    train_pred= model.predict_proba(train[cols])[:,1]

    density = True

    fig, ax = plt.subplots(figsize=(10, 7))

    background_color = 'red'

    opts = dict(
        range=[0,1],
        bins = 25,
        density = density
    )
    histtype1 = dict(
        histtype='stepfilled',
        linewidth=3,
        alpha=0.45,
    )

    ax.hist(train_background, **opts, **histtype1,
             facecolor=background_color,
             edgecolor=background_color,
             zorder=0)
    ax.hist(train_signal, **opts, **histtype1,
             facecolor='blue',
             edgecolor='blue',
             zorder=1000)

    hist_test_0 = np.histogram(test_background, **opts)
    hist_test_1 = np.histogram(test_signal, **opts)
    bins_mean = (hist_test_0[1][1:]+hist_test_0[1][:-1])/2
    bin_width = bins_mean[1]-bins_mean[0]
    area0 = bin_width*np.sum(test.label==0)
    area1 = bin_width*np.sum(test.label==1)

    opts2 = dict(
          capsize=3,
          ls='none',
          marker='o'
    )

    ax.errorbar(bins_mean, hist_test_0[0],  yerr = np.sqrt(hist_test_0[0]/area0), xerr=bin_width/2,
                 color=background_color, **opts2, zorder=100)
    ax.errorbar(bins_mean, hist_test_1[0],  yerr = np.sqrt(hist_test_1[0]/area1), xerr=bin_width/2,
                 color='blue', **opts2, zorder=10000)

    _ks_back = ks_2samp(train_background, test_background)[1]
    _ks_sign = ks_2samp(train_signal, test_signal)[1]



    auc_test  = roc_auc_score(test.label,test_pred )
    auc_train = roc_auc_score(train.label,train_pred)
    legend_elements = [Patch(facecolor='black', edgecolor='black', alpha=0.4,
                             label=f'Train (auc) : {round(auc_train,8)}'),
                      Line2D([0], [0], marker='|', color='black',
                             label=f'Test (auc) : {round(auc_test,8)}',
                              markersize=25, linewidth=1),
                       Circle((0.5, 0.5), radius=2, color='red',
                              label=f'Background (ks-pval) : {round(_ks_back,8)}',),
                       Circle((0.5, 0.5), 0.01, color='blue',
                              label=f'Signal (ks-pval) : {round(_ks_sign,8)}',),
                       ]

    ax.legend(
              #title='KS test',
              handles=legend_elements,
              #bbox_to_anchor=(0., 1.02, 1., .102),
              loc='upper center',
              ncol=2,
              #mode="expand",
              #borderaxespad=0.,
              frameon=True,
              fontsize=15)

    if print_params:
      ax.text(1.02, 1.02, params_to_string(model),
        transform=ax.transAxes,
      fontsize=13, ha='left', va='top')

    ax.set_yscale('log')
    ax.set_xlabel('XGB output')
    ax.set_ylim(0.005, 100)

    #plt.savefig(os.path.join(dir_, 'LR_overtrain.pdf'), bbox_inches='tight')
    return fig, ax

def roc(test_x,test_y,train_x,train_y, model):
    """"
    Presenta la curva roc, que muestra la precisión del clasificador, entre mas
    cercana sea el area 1 mejor será el clasificador """
    plt.figure(figsize=(10,7))
    plt.title('ROC curve', fontsize=20)
    model_predict = model.predict_proba(test_x)
    model_predict = model_predict[:,1]
    auc_score = roc_auc_score(test_y, model_predict)
    fpr, tpr, _ = roc_curve(test_y, model_predict) #roc_curve(true binary labels, prediction scores)
    print('Test : ', auc_score)
    plt.plot(tpr, 1-fpr, label='Test   '+ str(round(auc_score, 4)), color='purple', linewidth=3)

    model_predict = model.predict_proba(train_x)
    model_predict = model_predict[:,1]
    auc_score = roc_auc_score(train_y, model_predict)
    fpr, tpr, _ = roc_curve(train_y, model_predict)
    plt.plot(tpr, 1-fpr, label='Train   ' + str(round(auc_score,4)) , color='orange', linewidth=3)
    print('Train : ', auc_score)
    plt.legend(loc='best',fontsize=20)
    plt.ylabel('Background Rejection', fontsize=20)
    plt.xlabel('Signal efficiency', fontsize=20)
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.ylim(0.7,1)
    plt.show()



def stable_cumsum(arr, axis=None, rtol=1e-05, atol=1e-08):
    """Use high precision for cumsum and check that final value matches sum.

    Warns if the final cumulative sum does not match the sum (up to the chosen
    tolerance).

    Parameters
    ----------
    arr : array-like
        To be cumulatively summed as flat.
    axis : int, default=None
        Axis along which the cumulative sum is computed.
        The default (None) is to compute the cumsum over the flattened array.
    rtol : float, default=1e-05
        Relative tolerance, see ``np.allclose``.
    atol : float, default=1e-08
        Absolute tolerance, see ``np.allclose``.

    Returns
    -------
    out : ndarray
        Array with the cumulative sums along the chosen axis.
    """
    out = np.cumsum(arr, axis=axis, dtype=np.float64)
    expected = np.sum(arr, axis=axis, dtype=np.float64)
    if not np.allclose(
        out.take(-1, axis=axis), expected, rtol=rtol, atol=atol, equal_nan=True
    ):
        warnings.warn(
            (
                "cumsum was found to be unstable: "
                "its last element does not correspond to sum"
            ),
            RuntimeWarning,
        )
    return out


def get_true_false_positives(y_true, y_score):
    """
    Return True and False positive counts for each threshold value.

    Parameters:
    - y_true: array-like, true labels (0 or 1)
    - y_score: array-like, predicted scores or probabilities

    Returns:
    - true positives: array-like, 
    - false positives: array-like,
    - thresholds: array-like, threshold values
    """
    tps = []
    fns = []
    y_true_ = np.array(y_true)
    y_score_ = np.array(y_score)

    desc_score_indices = np.argsort(y_score_, kind="mergesort")[::-1]
    y_score_ = y_score_[desc_score_indices]
    y_true_ = y_true_[desc_score_indices]
    # y_score typically has many tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate a value for the end of the curve.
    distinct_value_indices = np.where(np.diff(y_score_))[0]
    threshold_idxs = np.r_[distinct_value_indices, y_true_.size - 1]

    # accumulate the true positives with decreasing threshold
    tps = stable_cumsum(y_true_)[threshold_idxs]
    #fps = 1 + threshold_idxs - tps
    fns = stable_cumsum(1-y_true_)[threshold_idxs]

    return tps, fns, y_score_[threshold_idxs]


def fast_max_fom(y_true, y_score):
    tps, fns, thr = get_true_false_positives(y_true, y_score)
    eff = tps/np.sum(y_true)
    fom = eff/np.sqrt(fns)

    mask = np.isfinite(fom)
    max_fom = np.max(fom[mask])
    threshold = thr[mask][np.argmax(fom[mask])]

    return max_fom, threshold


def correlation_heatmap(df, **kwargs):
    correlations = df.corr()
    heatmap(correlations, **kwargs)
    # fig, ax = plt.subplots(figsize=(10,10))
    # heatmap_args = dict( vmax=1.0, vmin=-1, center=0, fmt='.2f', 
    #                      #cmap="YlGnBu", 
    #                      cmap="coolwarm", 
    #                      square=True, linewidths=.5, cbar_kws={"shrink": .70}, 
    #                      annot=True, annot_kws={"size": 35 / np.sqrt(len(correlations))})
    # heatmap_args.update(kwargs)
    # sns.heatmap(correlations,
    #              **heatmap_args 
    #             )
    # plt.xticks(rotation=45)

def heatmap(df, **kwargs):
    heatmap_args = dict( vmax=1.0, vmin=-1, center=0, fmt='.2f', 
                         #cmap="YlGnBu", 
                         cmap="coolwarm", 
                         square=True, linewidths=.5, cbar_kws={"shrink": .90}, 
                         annot=True, annot_kws={"size": 35 / np.sqrt(len(df))})
    heatmap_args.update(kwargs)
    sns.heatmap(df,**heatmap_args )
    plt.xticks(rotation=45)

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def open_json(path):
    if not path:
        return dict()
    with open(path, 'r') as jj:
        return json.load(jj)


def save_json(path, obj):
    if not path: raise NotImplementedError
    with open(path, 'w') as jj:
        return json.dump(obj, jj, indent=4, cls=NpEncoder)

def extract_highly_correlated_pairs(df, threshold):
    """
    Extract feature pairs with correlation greater than the specified threshold.

    Parameters:
    - df: Pandas DataFrame
        The input DataFrame containing numerical features.
    - threshold: float
        The correlation threshold.

    Returns:
    - List of Lists
        A list of lists containing feature pairs and their correlation values.
    """
    # Calculate the correlation matrix
    corr_matrix = df.corr()

    # Extract feature pairs with correlation greater than the threshold
    correlated_pairs = []
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if abs(np.abs(corr_matrix.iloc[i, j])) > threshold:
                feature_A = corr_matrix.columns[i]
                feature_B = corr_matrix.columns[j]
                correlation_value = corr_matrix.iloc[i, j]
                correlated_pairs.append([feature_A, feature_B, correlation_value])

    return correlated_pairs



def eval_cross_validations(model, train_df, folds, vars=list(),
                            auc_=True, f1_50=True, avg_precision=True,
                            fast_fom=True,
                            verbose=False):
    """
    Evaluate a model using cross-validation and calculate specified metrics.

    Parameters:
    - model: The machine learning model to evaluate.
    - train_df: The training DataFrame containing features and labels.
    - folds: A list of tuples representing the train-test splits for cross-validation.
    - vars: A list of feature variables to use in the evaluation.
    - auc_: If True, calculate AUC score.
    - f1_50: If True, calculate F1 score with a threshold of 0.5.
    - avg_precision: If True, calculate average precision score.
    - verbose: If True, print intermediate results after each fold.

    Returns:
    - dict_var: A dictionary containing evaluation metrics for each fold.
      Format: {'metric': {'train': [values], 'test': [values]}, ...}
    """
    metrics = list()
    if auc_: metrics.append('auc')
    if f1_50: metrics.append('f1_50')
    if avg_precision: metrics.append('avg_precision')
    if fast_fom: metrics.append('FOM')
    
    dict_var = dict()
    
    # Initialize the dictionary structure for each metric
    for metric in metrics:
        dict_var[metric] = dict(train=list(), test=list())

    # Loop through the folds for cross-validation
    for fold_train, fold_test in folds:
        # Extract the training and testing subsets
        fold_train_df = train_df.iloc[fold_train]
        fold_test_df  = train_df.iloc[fold_test]

        evaluation = [(fold_train_df[vars], fold_train_df.label),
                        (fold_test_df[vars],  fold_test_df.label)]
        
        # Fit the model on the training data
        model.fit(fold_train_df[vars], 
                    fold_train_df.label,
                    eval_set=evaluation,
                    verbose=False,)
        
        # Predict probabilities for training and testing data
        proba_train = model.predict_proba(fold_train_df[vars])[:,1]
        proba_test  = model.predict_proba(fold_test_df[vars])[:,1]
        
        # Calculate and store metrics for each specified evaluation metric
        if auc_:
            train_auc = roc_auc_score(fold_train_df.label, proba_train)
            test_auc  = roc_auc_score(fold_test_df.label, proba_test)
            dict_var['auc']['train'].append(train_auc)
            dict_var['auc']['test'].append(test_auc) 
        
        if f1_50:
            train_f1 = f1_score(fold_train_df.label, proba_train>0.5, average='weighted')
            test_f1  = f1_score(fold_test_df.label, proba_test>0.5, average='weighted')
            dict_var['f1_50']['train'].append(train_f1)
            dict_var['f1_50']['test'].append(test_f1) 

        if avg_precision:
            train_avg = average_precision_score(fold_train_df.label, proba_train)
            test_avg  = average_precision_score(fold_test_df.label, proba_test)
            dict_var['avg_precision']['train'].append(train_avg)
            dict_var['avg_precision']['test'].append(test_avg) 

        if fast_fom:
            train_fom = fast_max_fom(fold_train_df.label, proba_train)
            test_fom  = fast_max_fom(fold_test_df.label, proba_test)
            dict_var['FOM']['train'].append(train_fom[0])
            dict_var['FOM']['test'].append(test_fom[0]) 
            

        # Print intermediate results if verbose is True
    if verbose:
        pprint(dict_var)

    return dict_var