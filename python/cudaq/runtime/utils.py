# ============================================================================ #
# Copyright (c) 2022 - 2024 NVIDIA Corporation & Affiliates.                   #
# All rights reserved.                                                         #
#                                                                              #
# This source code and the accompanying materials are made available under     #
# the terms of the Apache License 2.0 which accompanies this distribution.     #
# ============================================================================ #
from __future__ import annotations

from ..mlir._mlir_libs._quakeDialects import cudaq_runtime
from ..kernel.kernel_builder import PyKernel
from ..kernel.kernel_decorator import PyKernelDecorator
from ..mlir.dialects import quake, cc

import numpy as np
import sys
from typing import List


def __getArgTypes(kernel):
    argTypes = None
    if isinstance(kernel, PyKernel):
        argTypes = kernel.mlirArgTypes
    elif isinstance(kernel, PyKernelDecorator):
        argTypes = kernel.signature

    return argTypes


def __isFirstArgTypeStdVec(kernel, argTypes):
    isFirstArgTypeStdVec = None

    if isinstance(kernel, PyKernel):
        isFirstArgTypeStdVec = cc.StdvecType.isinstance(argTypes[0])
    elif isinstance(kernel, PyKernelDecorator):
        firstArgType = next(iter(argTypes))
        checkList = [
            list,
            np.ndarray,
            List,
            List[float],
            List[complex],
            List[int],
            'list',
            'np.ndarray',
            'List',
            'List[float]',
            'List[complex]',
            'List[int]',
            ## [PYTHON_VERSION_FIX] sys.version_info >= (3, 9)
            list[float],
            list[complex],
            list[int],
            list[bool],
            'list[float]',
            'list[complex]',
            'list[int]',
            'list[bool]'
        ]
        isFirstArgTypeStdVec = argTypes[firstArgType] in checkList

    return isFirstArgTypeStdVec


def __isBroadcast(kernel, *args):
    # kernel could be a PyKernel or PyKernelDecorator
    argTypes = __getArgTypes(kernel)
    if not argTypes or not args:
        return False

    # Quick check, if we have a 2d array anywhere, we know this is a broadcast
    isDefinitelyBroadcast = any(
        hasattr(arg, "shape") and len(arg.shape) == 2 for arg in args)

    if isDefinitelyBroadcast:
        # Error check, did the user pass a single value for any of the other arguments
        for i, arg in enumerate(args):
            if isinstance(arg, (int, float, bool, str)):
                raise RuntimeError(
                    f"2D array argument provided for an sample or observe broadcast, but argument {i} ({type(arg)}) must be a list."
                )

    firstArg = args[0]
    firstArgTypeIsStdvec = __isFirstArgTypeStdVec(kernel, argTypes)

    if isinstance(firstArg, (list, List)) and not firstArgTypeIsStdvec:
        return True

    if hasattr(firstArg, "shape"):
        if firstArg.ndim == 1 and not firstArgTypeIsStdvec:
            return True

        if firstArg.ndim == 2:
            return True

    return False


def __createArgumentSet(*args):
    nArgSets = len(args[0])
    argSet = [
        tuple(arg[j].tolist() if hasattr(arg, "tolist") and arg.ndim == 2 else
              arg[j] if isinstance(arg, (list, List)) else arg[j]
              for arg in args)
        for j in range(nArgSets)
    ]
    return argSet
