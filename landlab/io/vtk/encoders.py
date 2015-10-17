#!/bin/env python

import base64
import numpy as np


class EncoderError(Exception):
    pass


class UnknownEncoderError(EncoderError):

    def __init__(self, name):
        self._name

    def __str__(self):
        return '%s: Unknown encoder' % self._name


class EncoderInterface(object):

    def encode(self, array):
        pass


class AsciiEncoder(object):

    def encode(self, array):
        try:
            return ' '.join([str(val) for val in array.flatten()])
        except AttributeError:
            return ' '.join([str(val) for val in array])


class RawEncoder(object):

    def encode(self, array):
        try:
            as_str = array.tostring()
        except AttributeError:
            as_str = np.array(array).tostring()
        block_size = np.array(len(as_str), dtype=np.int32).tostring()
        return block_size + as_str


class Base64Encoder(object):

    def encode(self, array):
        try:
            as_str = array.tostring()
        except AttributeError:
            as_str = np.array(array).tostring()
        as_str = base64.b64encode(as_str)
        block_size = base64.b64encode(np.array(len(as_str),
                                               dtype=np.int32).tostring())
        return block_size + as_str

    def decode(self, array):
        pass


_ENCODERS = {
    'ascii': AsciiEncoder(),
    'raw': RawEncoder(), 'base64': Base64Encoder()
}


def encode(array, encoding='ascii'):
    return _ENCODERS[encoding].encode(array)


def decode(array, encoding='raw'):
    return _ENCODERS[encoding].decode(array)
