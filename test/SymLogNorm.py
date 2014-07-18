"""
This was taken from the source code of Matplotlib because it
is a new feature of the library and this helps with backwards
compatibility.
"""

from matplotlib.colors import *

class SymLogNorm(Normalize):
    """
    The symmetrical logarithmic scale is logarithmic in both the
    positive and negative directions from the origin.

    Since the values close to zero tend toward infinity, there is a
    need to have a range around zero that is linear.  The parameter
    *linthresh* allows the user to specify the size of this range
    (-*linthresh*, *linthresh*).
    """
    def __init__(self,  linthresh, linscale=1.0,
                 vmin=None, vmax=None, clip=False):
        """
        *linthresh*:
        The range within which the plot is linear (to
        avoid having the plot go to infinity around zero).

        *linscale*:
        This allows the linear range (-*linthresh* to *linthresh*)
        to be stretched relative to the logarithmic range.  Its
        value is the number of decades to use for each half of the
        linear range.  For example, when *linscale* == 1.0 (the
        default), the space used for the positive and negative
        halves of the linear range will be equal to one decade in
        the logarithmic range. Defaults to 1.
        """
        Normalize.__init__(self, vmin, vmax, clip)
        self.linthresh = float(linthresh)
        self._linscale_adj = (linscale / (1.0 - np.e ** -1))

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)
        vmin, vmax = self.vmin, self.vmax

        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            result.fill(0)
        else:
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)
            # in-place equivalent of above can be much faster
            resdat = self._transform(result.data)
            resdat -= self._lower
            resdat /= (self._upper - self._lower)

        if is_scalar:
            result = result[0]
        return result

    def _transform(self, a):
        """
        Inplace transformation.
        """
        masked = np.abs(a) > self.linthresh
        sign = np.sign(a[masked])
        log = (self._linscale_adj + np.log(np.abs(a[masked]) / self.linthresh))
        log *= sign * self.linthresh
        a[masked] = log
        a[~masked] *= self._linscale_adj
        return a

    def _inv_transform(self, a):
        """
        Inverse inplace Transformation.
        """
        masked = np.abs(a) > (self.linthresh * self._linscale_adj)
        sign = np.sign(a[masked])
        exp = np.exp(sign * a[masked] / self.linthresh - self._linscale_adj)
        exp *= sign * self.linthresh
        a[masked] = exp
        a[~masked] /= self._linscale_adj
        return a

    def _transform_vmin_vmax(self):
        """
        Calculates vmin and vmax in the transformed system.
        """
        vmin, vmax = self.vmin, self.vmax
        arr = np.array([vmax, vmin]).astype(np.float)
        self._upper, self._lower = self._transform(arr)

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        val = ma.asarray(value)
        val = val * (self._upper - self._lower) + self._lower
        return self._inv_transform(val)

    def autoscale(self, A):
        """
        Set *vmin*, *vmax* to min, max of *A*.
        """
        self.vmin = ma.min(A)
        self.vmax = ma.max(A)
        self._transform_vmin_vmax()

    def autoscale_None(self, A):
        """ autoscale only None-valued vmin or vmax """
        if self.vmin is not None and self.vmax is not None:
            pass
        if self.vmin is None:
            self.vmin = ma.min(A)
        if self.vmax is None:
            self.vmax = ma.max(A)
        self._transform_vmin_vmax()