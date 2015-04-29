import os
from sys import exit

from .log import get_logger
logger = get_logger()

def get_workspace():
    from ROOT import RooWorkspace
    return RooWorkspace('workspace')

def filter_content(lines):
    lineno = 0
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            yield (line, lineno)
        lineno += 1

def assemble_model(file):
    """
    :param file: a filename for a .model file or an iterable
    :returns: a workspace filled with the model described in `file`
    """
    if type(file) == str:
        with open(file) as f:
            return assemble_model(f)

    w = get_workspace()
    factory = w.factory()

    logger.info('Assembling model file: \'%s\'' % file)

    for line, lineno in filter_content(file):
        if not factory.process(line):
            logger.error('RooFit factory error at line %d', lineno)
            exit(1)

    return w

def split_model(model, dataset, split_on, differing):
    """
    :param model:
    :param dataset:
    :param splitOn:
    :param model_name:
    :param differing:
    :returns: A simultaneous pdf split over the provided category
    """
    from ROOT import RooArgSet, RooSimPdfBuilder, SetOwnership

    model_name = model.GetName()

    # physModels = model1 [model2..]
    # splitCats  = splitVar1 [splitVar2..]
    # model1     = splitVar1 : diff1 [diff2..]

    logger.debug('split_model parameters are %s, %s and %s' % (split_on, model_name, differing))

    builder = RooSimPdfBuilder(RooArgSet(model))
    bconfig = builder.createProtoBuildConfig()
    bconfig['physModels'] =  model_name
    bconfig['splitCats'] = split_on
    bconfig[model_name] = "%s : %s" % (split_on, ','.join(differing))

    pdf = builder.buildPdf(bconfig, dataset)

    SetOwnership(builder, False)

    return pdf

def split_custom(models, dataset, config):
    """
    :config: Dictionary with the following format:
             { 'model1': [('splitCat1', ['var1', 'var2', 'var3']),
                          ('splitCat2', ['var2', 'var3'        ]),
                         ],
               'model2': [('splitCat1', ['var3'])]
             }
    """

    from ROOT import RooArgSet, RooSimPdfBuilder, SetOwnership
    model_name = model.GetName()
    logger.debug('split_custom configuration is %s' % config)

    builder = RooSimPdfBuilder(RooArgSet(*models))
    bconfig = builder.createProtoBuildConfig()

    models = config.keys()
    splitCats = []
    for m in config:
        splitCats.append(config[m][0])

    bconfig['physModels'] = ','.join(config.keys())
    bconfig['splitCats'] = ','.join(config.keys())

    for m in config:
        cats = bconfig[m]
        entries = ""
        for c in cats:
            entries.append('%s : %s\n' % (c[0], ','.join(c[1])))

        bconfig[m] = '\n'.join(entries)

    pdf = builder.buildPdf(bconfig, dataset)
    SetOwnership(builder, False)

    return pdf

def get_values(model, variables):
    from ROOT import RooArgSet
    result = [None] * len(variables)
    result_err = [None] * len(variables)
    it = model.getParameters(RooArgSet()).iterator()
    next = it.Next()
    while next:
        if next.GetName() in variables:
            idx = variables.index(next.GetName())
            result[idx] = next.getVal()
            result_err[idx] = next.getError()
        next = it.Next()
    if None in variables:
        logger.error('Variable not found in model: %s' % (variables[result.index(None)]))
        exit(1)
    return result, result_err

def mle(model, data, start_params='None', out_params='result.params', numcpus=1, strategy=1, fit=True, extended=False):
    from ROOT import RooFit
    if start_params == 'None':
        logger.warning('Not reading in initial parameters. Using defaults from model.')
    elif not os.path.exists(start_params):
        logger.error('No such parameter file: \'%s\'' % start_params)
        exit(1)
    else:
        logger.info('Reading in parameters from \'%s\'' % start_params)
        model.getParameters(data).readFromFile(start_params)
    logger.info('Starting fit')
    if fit:
        results = model.fitTo(
                data,
                RooFit.NumCPU(numcpus),
                RooFit.Extended(extended),
                RooFit.Minos(False),
                RooFit.Strategy(2),
                RooFit.Minimizer('Minuit2'),
                RooFit.Save(True),
        )
        logger.info('Finished fit')
    else:
        logger.info('Not performing fit')
    logger.info('Writing resulting parameters to \'%s\'' % out_params)
    model.getParameters(data).writeToFile(out_params)
    return results

def add_weights(model, data, yield_names):
    iter = model.getParameters(data).iterator()
    yields = []
    next = iter.Next()
    while next:
        if not next.GetName() in yield_names:
            next.setConstant(True)
        else:
            yields.append(next)
        next = iter.Next()

    from ROOT import RooStats, RooDataSet, RooArgList
    sdata = RooStats.SPlot("SPlot", "SPlot", data, model, RooArgList(*yields))
    components = []
    for name in yield_names:
        components.append(RooDataSet(data.GetName(), data.GetTitle(), data, data.get(), '', name + '_sw'))
    return components

def get_matching_vars(wspace, tree):
    varSet = wspace.allVars()
    vars = dict()
    it = varSet.createIterator()
    el = it.Next()
    result = []
    while el:
        vars[el.GetName()] = el
        el = it.Next()
    catSet = wspace.allCats()
    it = catSet.createIterator()
    el = it.Next()
    while el:
        vars[el.GetName()] = el
        el = it.Next()

    for b in tree.GetListOfBranches():
        if b.GetName() in vars:
            result.append(vars[b.GetName()])

    logger.info('Matched the following variables between model and data: %s' % [v.GetName() for v in result])
    return result

def load_tree(w, fname, treename, cutstring=''):
    from ROOT import RooDataSet, TFile, RooArgSet

    if not os.path.exists(fname):
        raise IOError('File not found: %s' % fname)

    logger.debug('Loading data from tree: %s/%s' % (fname, treename))

    file = TFile(fname)
    tree = file.Get(treename)

    if cutstring:
        tmp = TFile('tmp.root', 'RECREATE')
        logger.info('Applying cut string: \'%s\'' % cutstring)
        tree = tree.CopyTree(cutstring)

    vars = get_matching_vars(w, tree)
    tree.SetBranchStatus("*", False)
    for v in vars:
        tree.SetBranchStatus(v.GetName(), True)
    args = RooArgSet(*vars)
    data = RooDataSet('data', 'data', tree, args)

    return data

def add_weights(model, data, yield_names):
    iter = model.getParameters(data).iterator()
    yields = []
    next = iter.Next()
    while next:
        if not next.GetName() in yield_names:
            next.setConstant(True)
        else:
            yields.append(next)
        next = iter.Next()

    from ROOT import RooStats, RooDataSet, RooArgList
    sdata = RooStats.SPlot("SPlot", "SPlot", data, model, RooArgList(*yields))
    #components = []
    #for name in yield_names:
    #    components.append(RooDataSet(data.GetName(), data.GetTitle(), data, data.get(), '', name + '_sw'))
    #return components

