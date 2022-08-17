#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is the script to convert old json files to new one.
"""

__author__ = "Gene Ting-Chun Kao"
__email__ = "kao@arch.ethz.ch"


if __name__ == '__main__':
    import compas
    import compas_cra
    import os

    from compas_cra.datastructures import CRA_Assembly

    file_name = './type-d.json'
    assembly = CRA_Assembly.from_json(
        os.path.join(compas_cra.DATA, file_name))

    compas.json_dump(assembly,
                     os.path.join(compas_cra.DATA, file_name))

    print("converted")
