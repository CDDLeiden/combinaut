from scaffviz.clustering.manifold import Manifold


class Manifold(Manifold):

    def __init__(self, method):
        self.method = method

    def fit(self, X):
        self.method.fit(X)
        return self

    def transform(self, X):
        return self.method.fit_transform(X)

    def fit_transform(self, X):
        return self.method.fit_transform(X)

    def __str__(self):
        return f"{self.method.__class__.__name__}"
