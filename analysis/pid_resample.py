
class Resampler:
    def __init__(self, *args):
        # Choose histogram size according to bin edges
        # Take under/overflow into account for dependent variables only
        edges = []
        for arg in args[:-1]:
            edges.append(np.append(np.append([-np.inf], arg), [np.inf]))
        edges.append(args[-1])
        self.edges = edges

        self.histogram = np.zeros(map(lambda x: len(x) - 1, self.edges))

    def learn(self, features, weights=None):
        assert(len(features) == len(self.edges))

        features = np.array(features)

        h , _ = np.histogramdd(features.T, bins=self.edges, weights=weights)
        self.histogram += h

    def sample(self, features):
        import numpy as np

        assert(len(features) == len(self.edges) - 1)

        args = np.array(features)

        idx = [np.searchsorted(edges, vals) - 1 for edges, vals in zip(self.edges, args)]

        tmp = self.histogram[idx]
        # Fix negative bins (resulting from possible negative weights) to zero
        tmp[tmp < 0] = 0

        norm = np.sum(tmp, axis=1)
        probs = tmp / norm[:,np.newaxis]

        sampled_bin = []
        for i in range(tmp.shape[0]):
            sampled_bin.append(np.random.choice(tmp.shape[1], p=probs[i,:]))
        sampled_bin = np.array(sampled_bin)
        sampled_val = np.random.uniform(self.edges[-1][sampled_bin],
                                        self.edges[-1][sampled_bin + 1],
                                        size=len(sampled_bin))

        # If the histogram is empty, we can't sample
        sampled_val[norm == 0] = np.nan

        return sampled_val
