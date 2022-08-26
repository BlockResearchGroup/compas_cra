from __future__ import print_function

import platform
import compas_cra


if __name__ == "__main__":

    print()
    print("Yay! COMPAS CRA is installed correctly!")
    print()
    print("COMPAS CRA: {}".format(compas_cra.__version__))
    print(
        "Python: {} ({})".format(
            platform.python_version(), platform.python_implementation()
        )
    )
