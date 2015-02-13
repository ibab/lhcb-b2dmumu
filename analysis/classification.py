from joblib import Parallel, delayed

import logging

N_trees = 1000

def run_crossval(clf, i, X_train, y_train, X_test, y_test):
    import numpy as np
    from sklearn.metrics import roc_curve, zero_one_loss
    logging.info('    Running fold #{}'.format(i + 1))
    probs = clf.fit(X_train, y_train).predict_proba(X_test)
    fpr, tpr, thresholds = roc_curve(y_test, probs[:,1])

    err = np.zeros((N_trees,))
    for j, y_pred in enumerate(clf.staged_predict(X_test)):
        err[j] = zero_one_loss(y_pred, y_test)

    logging.info('    Finished fold #{}'.format(i + 1))

    return fpr, tpr, err


def validate_classifier(clf, X, y, outputdir):
    import numpy as np
    from sklearn.cross_validation import StratifiedKFold
    from sklearn.metrics import roc_curve, auc, zero_one_loss
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    sns.set_palette('deep', desat=.6)
    sns.set_context('talk')

    with PdfPages(outputdir + '/classifier.pdf') as pdf:
        logging.info('Calculating correlation...')
        import pandas
        #corr = pandas.DataFrame(sidebands).corr()
        #sns.heatmap(corr, vmax=.8, linewidths=0, square=True)
        #plt.tight_layout()
        #pdf.savefig()
        #plt.clf()

        logging.info('Running x-validation...')
        skf = StratifiedKFold(y, 5)

        mean_fpr = np.linspace(0, 1, 200)
        mean_tpr = np.zeros(200)

        errs = []

        results = Parallel(n_jobs=5)(delayed(run_crossval)(clf, i, X[train], y[train], X[test], y[test]) for i, (train, test) in enumerate(skf))
        for f, t, e in results:
            plt.plot(f, t, lw=1)

        #for i, (train, test) in enumerate(skf):
        #    logging.info('    Running fold #{}'.format(i + 1))
        #    probs = clf.fit(X[train], y[train]).predict_proba(X[test])
        #    fpr, tpr, thresholds = roc_curve(y[test], probs[:,1])
        #    plt.plot(fpr, tpr, lw=1)

        #    err = np.zeros((N_trees,))
        #    for i, y_pred in enumerate(clf.staged_predict(X[test])):
        #        err[i] = zero_one_loss(y_pred, y[test])
        #    errs.append(err)

        plt.ylim(0.8, 1.0)
        plt.xlim(0.0, 0.2)
        plt.ylabel('True positive rate')
        plt.xlabel('False positive rate')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        for f, t, e in results:
            plt.plot(np.arange(N_trees)+1, e)
        plt.ylabel('Error')
        plt.xlabel('$N_\\mathrm{estimators}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        #logging.info('Plot importances...')
        #imp = sorted(zip(sidebands.dtype.names, clf.feature_importances_), key=lambda x: x[1])
        #plt.barh(np.arange(len(imp)), [entr[1] for entr in imp], 0.30, align='center')
        #plt.yticks(np.arange(len(imp)), [entr[0] for entr in imp])
        #plt.tight_layout()
        #pdf.savefig()
        #plt.clf()

