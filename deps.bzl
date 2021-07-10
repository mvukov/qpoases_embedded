""" Configures repo dependencies.
"""

load("@bazel_skylib//:workspace.bzl", "bazel_skylib_workspace")

def qpoases_embedded_deps():
    bazel_skylib_workspace()
