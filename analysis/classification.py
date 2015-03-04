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

    err_train = np.zeros((N_trees,))
    for j, y_pred in enumerate(clf.staged_predict(X_train)):
        err_train[j] = zero_one_loss(y_pred, y_train)

    logging.info('    Finished fold #{}'.format(i + 1))

    return fpr, tpr, err, err_train


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

        plt.figure(figsize=(10, 10))
        results = Parallel(n_jobs=5)(delayed(run_crossval)(clf, i, X[train], y[train], X[test], y[test]) for i, (train, test) in enumerate(skf))
        for f, t, _, _ in results:
            plt.plot(1 - f, t, lw=2)

        #for i, (train, test) in enumerate(skf):
        #    logging.info('    Running fold #{}'.format(i + 1))
        #    probs = clf.fit(X[train], y[train]).predict_proba(X[test])
        #    fpr, tpr, thresholds = roc_curve(y[test], probs[:,1])
        #    plt.plot(fpr, tpr, lw=1)

        #    err = np.zeros((N_trees,))
        #    for i, y_pred in enumerate(clf.staged_predict(X[test])):
        #        err[i] = zero_one_loss(y_pred, y[test])
        #    errs.append(err)

        plt.ylim(0.6, 1.0)
        plt.xlim(0.8, 1.0)
        plt.ylabel('True positive rate')
        plt.xlabel('True negative rate')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        import seaborn as sns
        sns.set_context('talk')
        for f, t, e, e_train in results:
            plt.plot(np.arange(N_trees)+1, e)
            plt.plot(np.arange(N_trees)+1, e_train, alpha=0.5)
        plt.ylabel('Zero-one loss')
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

