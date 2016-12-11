class BaseVectorFunction(object):
    """
    This class represents a base vector function
    """
    def __init__(self, function_vector, **kwargs):
        """
        Initialize the vector function
        :param function_vector: - an iterable of two-argument lambda functions
        :return : - None
        """
        self.function_vector = function_vector


class DoubleVectorFunction(BaseVectorFunction):
    """
    This class represents a vector function of type F(t, U), namely a function
    of a single scalar and a single vector argument
    """

    def __call__(self, t, U):
        return [func(t, U) for func in self.function_vector]


class VectorFunction(BaseVectorFunction):
    """
    This class represents a vector function of type U(t), namely a vector
    function of a single scalar argument
    """

    def __call__(self, t):
        return [func(t) for func in self.function_vector]
