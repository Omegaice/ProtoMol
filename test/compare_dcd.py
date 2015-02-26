#!/usr/bin/python

import os
import math
import logging
from struct import unpack

def absDiff(one, two, ignoreSign=False):
    if ignoreSign:
        (one, two) = (math.fabs(one), math.fabs(two))
    return max(one, two) - min(one, two)

def compare_dcd(fExpected, fNew, epsilon, scalar=1.0, ignoreSign=False):
    dcdExpected = DCDReader(fExpected)
    dcdExpected.open()
    fExpected = dcdExpected.read_to_dict()
    dcdExpected.close()

    dcdActual = DCDReader(fNew)
    dcdActual.open()
    fActual = dcdActual.read_to_dict()
    dcdActual.close()

    if fExpected['cord'] != fActual['cord']:
        logging.debug("Cords Differ. Should be %s but is %s" % (fExpected['cord'], fActual['cord']))
        return False

    if fExpected['frames'] != fActual['frames']:
        logging.debug("Frame Count Differs. Should be %d but is %d" % (fExpected['frames'], fActual['frames']))
        return False

    if fExpected['firststep'] != fActual['firststep']:
        logging.debug("First Step Differs. Should be %d but is %d" % (fExpected['firststep'], fActual['firststep']))
        return False

    #if fExpected['comment'] != fActual['comment']:
    #    logging.debug("Comment Differs. Should be %s but is %s" % (fExpected['comment'], fActual['comment']))
    #    return False

    if fExpected['numatoms'] != fActual['numatoms']:
        logging.debug("Atom Count Differs. Should be %d but is %d" % (fExpected['numatoms'], fActual['numatoms']))
        return False

    diffs = 0
    for frame in xrange(0, len(fExpected['coordinates'])):
        for atom in xrange(0, fExpected['numatoms']):
            for element in xrange(0, 3):
                feone = float(fExpected['coordinates'][frame][atom][element]) * scalar
                fetwo = float(fActual['coordinates'][frame][atom][element])

                fediff = absDiff(feone, fetwo, ignoreSign)

                if fediff > epsilon:
                    diffs = diffs + 1
                    logging.debug('Frame %d, Atom %d, Element %d Differs' % (frame, atom, element))
                    logging.debug('Expected: %f, Actual: %f, Difference: %f' % (feone, fetwo, fediff))

    return diffs == 0

class DCDReader(object):
    def __init__(self, flname):
        self.flname = flname
        self.fl = None

    def open(self):
        self.fl = open(self.flname)

    def close(self):
        self.fl.close()
        self.fl = None

    def __read_int(self):
        buff = self.fl.read(4)
        return unpack("i", buff)[0]

    def __read_float(self):
        buff = self.fl.read(4)
        return unpack("f", buff)[0]

    def read_header(self):
        #endian_flag = self.__read_int()
        self.header_size = self.__read_int()
        self.cord = self.fl.read(4)
        self.frames = self.__read_int()
        self.firststep = self.__read_int()
        self.ignore1 = self.fl.read(24)
        self.freeIndexes = self.__read_int()
        self.ignore2 = self.fl.read(44)
        tail = self.__read_int()

        comment_length = self.__read_int()
        self.comment = self.fl.read(comment_length)
        comment_tail = self.__read_int()

        head = self.__read_int()
        self.numatoms = self.__read_int()
        tail = self.__read_int()

        if self.freeIndexes > 0:
            string = self.fl.read(4 * (self.numatoms - self.freeIndexes + 2))

        assert self.frames != 0, "Number of frames is 0!  Try fixing your DCD by running through catdcd"

    def __iter__(self):
        return self.read()

    def read(self):
        if self.fl == None:
            self.open()
        self.read_header()
        for i in xrange(self.frames):
            yield self.read_frame()
        self.close()

    def read_frame(self):
        head = self.__read_int()
        x = [self.__read_float() for i in xrange(self.numatoms)]
        tail = self.__read_float()
        head = self.__read_int()
        y = [self.__read_float() for i in xrange(self.numatoms)]
        tail = self.__read_int()
        head = self.__read_int()
        z = [self.__read_float() for i in xrange(self.numatoms)]
        tail = self.__read_int()
        return zip(x, y, z)

    def read_to_dict(self):
        self.read_header()
        d = dict(self.__dict__)
        del d["fl"]
        del d["flname"]
        d["coordinates"] = [self.read_frame() for i in xrange(self.frames)]
        return d
