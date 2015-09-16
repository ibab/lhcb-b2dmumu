
import numpy as np

def mksize(scale, aspect=(np.sqrt(5.0)-1)/2):
    import numpy as np
    fig_width_pt = 375.57643
    inches_per_pt = 1.0/72.27
    fig_width = fig_width_pt*inches_per_pt*scale
    fig_height = fig_width*aspect
    fig_size = [fig_width,fig_height]
    return fig_size


