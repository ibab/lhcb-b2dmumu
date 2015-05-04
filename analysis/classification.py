from joblib import Parallel, delayed

from analysis.log import get_logger
logger = get_logger()

N_trees = 1000

def run_crossval(clf, i, X_train, y_train, X_test, y_test, name, trees):
    import numpy as np
    from sklearn.metrics import roc_curve, zero_one_loss
    from sklearn.base import clone
    logger.info('    Running fold #{}'.format(i + 1))
    clf = clone(clf)
    logger.warn('after clone')
    probs = clf.fit(X_train, y_train).predict_proba(X_test)
    logger.warn('after fit')
    fpr, tpr, thresholds = roc_curve(y_test, probs[:,1])

    err = np.zeros((trees,))
    err_train = np.zeros((trees,))
    try:
        logger.info('    Calculating staged predict #{} (test)'.format(i + 1))
        for j, y_pred in enumerate(clf.staged_predict_proba(X_test)):
            p = y_pred[:,1]
            p[p > 0.5] = 1
            p[p <= 0.5] = 0
            err[j] = zero_one_loss(p, y_test)

        logger.info('    Calculating staged predict #{} (train)'.format(i + 1))
        for j, y_pred in enumerate(clf.staged_predict_proba(X_train)):
            p = y_pred[:,1]
            p[p > 0.5] = 1
            p[p <= 0.5] = 0
            err_train[j] = zero_one_loss(p, y_train)

    except AttributeError as e:
        logger.warn('staged_predict_proba not supported on {}'.format(name))

    sig_proba = clf.predict_proba(X_test[y_test == 1])[:,1]
    bkg_proba = clf.predict_proba(X_test[y_test == 0])[:,1]

    logger.info('    Finished fold #{}'.format(i + 1))
    results = {
            'fpr': fpr,
            'tpr': tpr,
            'err': err,
            'err_train': err_train,
            'sig_proba': sig_proba,
            'bkg_proba': bkg_proba,
            'y_test': y_test,
            'probs': probs[:,1],
            'clf': clf,
    }

    return results

def evaluate_classifier(clf, X, y, outputdir, name, folds=10, trees=1000):
    import numpy as np
    from sklearn.cross_validation import StratifiedKFold
    from sklearn.metrics import roc_curve, auc, zero_one_loss
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    #import seaborn as sns
    #sns.set_palette('deep', desat=.6)
    #sns.set_context('talk')

    with PdfPages(outputdir + '/' + name + '_classifier.pdf') as pdf:
        logger.info('Calculating correlation...')
        import pandas
        #corr = pandas.DataFrame(sidebands).corr()
        #sns.heatmap(corr, vmax=.8, linewidths=0, square=True)
        #plt.tight_layout()
        #pdf.savefig()
        #plt.clf()

        logger.info('Running x-validation...')
        skf = StratifiedKFold(y, folds)

        mean_fpr = np.linspace(0, 1, 200)
        mean_tpr = np.zeros(200)

        errs = []

        plt.figure(figsize=(10, 10))
        results = Parallel(n_jobs=folds)(delayed(run_crossval)(clf, i, X[train], y[train], X[test], y[test], name, trees) for i, (train, test) in enumerate(skf))

        total_y_test = []
        total_probs = []
        clfs = []
        for res in results:
            plt.plot(1 - res['fpr'], res['tpr'], lw=1, color='blue', alpha=0.6)
            total_y_test.append(res['y_test'])
            total_probs.append(res['probs'])
            clfs.append(res['clf'])

        fpr_, tpr_, thresh = roc_curve(np.hstack(total_y_test), np.hstack(total_probs))
        plt.plot(1 - fpr_, tpr_, lw=2, color='black')

        #for i, (train, test) in enumerate(skf):
        #    logger.info('    Running fold #{}'.format(i + 1))
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

        #import seaborn as sns
        #sns.set_context('talk')
        for res in results:
            plt.plot(np.arange(trees)+1, res['err'], color='blue', alpha=0.6)
            plt.plot(np.arange(trees)+1, res['err_train'], color='green', alpha=0.6)
        plt.ylabel('Zero-one loss')
        plt.xlabel('$N_\\mathrm{estimators}$')
        plt.tight_layout()
        pdf.savefig()
        plt.clf()

        total_signal = []
        total_bkg = []
        for res in results:
            total_signal.append(res['sig_proba'])
            total_bkg.append(res['bkg_proba'])

        def logit(x):
            return np.log(x / (1 - x))

        plt.hist(logit(np.hstack(total_signal)), range=(-20, 20), bins=100, alpha=0.6, color='blue', histtype='stepfilled', normed=True, label='Signal (Simulation)')
        plt.hist(logit(np.hstack(total_bkg)), range=(-20, 20), bins=100, alpha=0.6, color='red', histtype='stepfilled', normed=True, label='Background (Sidebands)')
        plt.xlabel('Classifier', ha='right', x=1)
        plt.ylabel('Candidates (normalized)', ha='right', y=1)
        plt.legend(loc='best')
        pdf.savefig()
        plt.clf()

        logger.warn(np.hstack(total_bkg))
        logger.warn('Assign from: {}'.format(len(np.hstack(total_bkg))))

        return logit(np.hstack(total_bkg)), logit(np.hstack(total_signal)), clfs

        #logger.info('Plot importances...')
        #imp = sorted(zip(sidebands.dtype.names, clf.feature_importances_), key=lambda x: x[1])
        #plt.barh(np.arange(len(imp)), [entr[1] for entr in imp], 0.30, align='center')
        #plt.yticks(np.arange(len(imp)), [entr[0] for entr in imp])
        #plt.tight_layout()
        #pdf.savefig()
        #plt.clf()

