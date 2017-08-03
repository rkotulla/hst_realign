#!/usr/bin/env python

import os
import sys

import numpy
import scipy.optimize

import pyds9
import argparse
import astLib.astWCS as astWCS
import scipy.spatial

from astropy.io import fits
from astroquery.vizier import Vizier
import astropy.coordinates
from astropy.coordinates import Angle, FK5

import astropy.units as u
from astropy.units import Quantity
from astroquery.gaia import Gaia
import astropy.wcs

def new_wcs(wcs, p):

    new_wcs = wcs.copy()

    new_wcs.header['CRPIX1'] = p[0]
    new_wcs.header['CRPIX2'] = p[1]

    cos_a = 1. #numpy.cos(p[2])
    sin_a = 0. #numpy.sin(p[2])

    new_wcs.header['CD1_1'] = wcs.header['CD1_1']*cos_a + wcs.header['CD1_2']*sin_a
    new_wcs.header['CD1_2'] = wcs.header['CD1_2']*cos_a - wcs.header['CD1_1']*sin_a
    new_wcs.header['CD2_1'] = wcs.header['CD2_1']*cos_a + wcs.header['CD2_2']*sin_a
    new_wcs.header['CD2_2'] = wcs.header['CD2_2']*cos_a - wcs.header['CD2_1']*sin_a
    new_wcs.updateFromHeader()

    return new_wcs

def wcs_fit(p, wcs, headers, xy, radec, cos_dec):
    # update the WCS

    for i, key in enumerate(headers):
        wcs.header[key] = p[i]
    wcs.updateFromHeader()
    ra_dec = numpy.array(wcs.pix2wcs(xy[:, 0], xy[:, 1]))

    # _wcs = new_wcs(wcs, p)
    # ra_dec = numpy.array(_wcs.pix2wcs(xy[:,0], xy[:,1]))
    # print "iteration: %d %d %f" % (p[0], p[1], numpy.degrees(p[2]))

    diff = (radec - ra_dec) * [cos_dec, 1.0]
    return diff.flatten()




if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Fix HST absolute alignment')

    parser.add_argument(
        'filename', type=str,
        metavar='input.fits',
        help='filename of input file')
    parser.add_argument(
        '-extname', type=str, metavar='extname',
        help='extension name', default='SCI',
    )
    args = parser.parse_args()

    hdulist= fits.open(args.filename)
    hdulist.info()

    headers = ['CRPIX1', 'CRPIX2',
               'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']

    #
    # open ds9
    #
    pyds9.ds9_xpans()
    ds9 = pyds9.DS9(target="hst_realign", start="-scale zscale -zoom 0.5", wait=25)


    ext = hdulist[args.extname]
    backup_hdr = ext.header.copy()
    wcs = astWCS.WCS(ext.header, mode="pyfits")
    cos_dec = numpy.cos(numpy.radians(ext.header['CRVAL2']))
    datafile = "tmp_%s.cat" % (args.extname)

    print "***************************************\n"*2
    print "Working on EXT %s\n" % (ext.name)
    print "***************************************\n"*2


    # center_ra, center_dec = wcs.getCentreWCSCoords()
    # print center_ra, center_dec
    # # display image in ds9
    # pass
    #
    # # save current extension as separate file
    tmp_file = "tmp_%s.fits" % (ext.name)
    fits.PrimaryHDU(data=ext.data, header=ext.header).writeto(tmp_file, clobber=True)

    ds9.set("frame delete all")

    ds9.set("frame 1")
    ds9.set("wcs skyformat degrees")
    ds9.set("file %s" % (tmp_file))
    #
    #
    #     ds9.set("frame 2")
    #     ds9.set("wcs skyformat degrees")
    #     cmd = "dsssao coord %f %f degrees size 45 45 arcmin" % (center_ra, center_dec)
    #     print cmd
    #     ds9.set(cmd)
    #     ds9.set("dsssao close")
    #
    #
    #     ds9.set("frame 1")
    #     ds9.set("lock frame wcs")
    ds9.set("catalog 2mass")
    ds9.set("catalog regions")
    ds9.set("catalog clear")
    ds9.set("catalog close")
    #
    # #
    # # Now get a set of point-pairs, the first indicating the position in
    # # the image, the second one in the WCS system
    # #
    all_pixel = []
    all_wcs = []
    #
    if (os.path.isfile(datafile)):
        merged = numpy.loadtxt(datafile)
        print "Loading %d position from previous run" % (merged.shape[0])
        for i in range(merged.shape[0]):
            all_pixel.append(merged[i, 0:2])
            all_wcs.append(merged[i, 2:4])
    #
    happy = False
    # if (not args.interactive):
    #     if (len(all_pixel) >= 1):
    #         np_all_pixel = numpy.array(all_pixel)
    #         np_all_wcs = numpy.array(all_wcs)
    #         p_init = [0] * len(headers)
    #         for i, key in enumerate(headers):
    #             p_init[i] = ext.header[key]
    #         fit = scipy.optimize.leastsq(wcs_fit,
    #                                      p_init,
    #                                      args=(wcs, headers, np_all_pixel,
    #                                            np_all_wcs, cos_dec),
    #                                      maxfev=1000,
    #                                      full_output=1)
    #         improved_wcs = fit[0]
    #         print "WCS correction: %d/%d pixels" % (
    #             improved_wcs[0] - p_init[0],
    #             improved_wcs[1] - p_init[1],
    #             #     numpy.degrees(improved_wcs[2]-p_init[2]),
    #         )
    #         for i, key in enumerate(headers):
    #             ext.header[key] = improved_wcs[i]
    #
    # else:
    while (not happy):

        while (True):

            break
            while (True):
                print "get pixel position"
                try:
                    _coords = ds9.get("imexam coordinate image")
                    # print _coords
                    coords_px = [float(_coords.split()[0]),
                                 float(_coords.split()[1])]
                    break
                except:
                    pass
                # print coords, coords_ali

            if (_coords == "0 0"):
                break


            while (True):
                print "select sky position"
                try:
                    _coords = ds9.get("imexam coordinate fk5")
                    coords_wcs = [float(_coords.split()[0]),
                                  float(_coords.split()[1])]
                    #   print _coords, coords_wcs
                    break
                except:
                    pass

            # if (coords == "0 0"):
            #     break

            # draw a little arrow from star to wcs position
            # for this to work, convert the ra/dec to x/y in the image frame
            wcs_xy = wcs.wcs2pix(coords_wcs[0], coords_wcs[1])
            # print coords_px
            # print wcs_xy

            ds9.set("regions system image")
            region_command = "regions command {vector(%f,%f,%f,%f) # color=blue}" % (
                    coords_px[0], coords_px[1],
                    wcs_xy[0]-coords_px[0],
                    wcs_xy[1]-coords_px[1],
            )
            region_command = "regions command {line(%f,%f,%f,%f) # color=blue}" % (
                    coords_px[0], coords_px[1],
                    wcs_xy[0], wcs_xy[1],
            )
            # print region_command
            ds9.set(region_command)

            all_pixel.append(coords_px)
            all_wcs.append(coords_wcs)

        check = "\n" #raw_input("ok")
        if (check == "r"):
            print "Resetting"
            all_pixel, all_wcs = [], []
            continue
        elif (check == "x"):
            sys.exit(0)
    #
    #
    #
    #
        #
        # Now we have a full set of x/y and ra/dec pairs, modify the CD and
        # CRPIX values to match the WCS
        #
        np_all_pixel = numpy.array(all_pixel)
        np_all_wcs = numpy.array(all_wcs)

        p_init = [0] * len(headers)
        for i,key in enumerate(headers):
            p_init[i] = ext.header[key]

        # p_init = [0] * len(headers)
        # p_init = [ext.header['CRPIX1'], ext.header['CRPIX2'], 0]
            #
            # for i,key in enumerate(headers):
            # p_init[i] = ext.header[key]

        print p_init

        fit = scipy.optimize.leastsq(wcs_fit,
                                 p_init,
                                 args=(wcs, headers, np_all_pixel, np_all_wcs, cos_dec),
                                 maxfev=1000,
                                 full_output=1)
        # print fit
        improved_wcs = fit[0]

        print "WCS correction: %d/%d pixels" % (
            improved_wcs[0]-p_init[0],
            improved_wcs[1]-p_init[1],
        #     numpy.degrees(improved_wcs[2]-p_init[2]),
        )

        #
        # Now update the header of the input image and save this CCD
        #
        # _wcs = new_wcs(wcs, improved_wcs)
        # for i,key in enumerate(headers):
        #     print ext.header[key], "-->", _wcs.header[key]
        #     ext.header[key] = _wcs.header[key]
        for i,key in enumerate(headers):
            # print ext.header[key], "-->", _wcs.header[key]
            ext.header[key] = improved_wcs[i]

        fixed_file = "tmp_%s.wcsfix.fits" % (ext.name)
        fits.PrimaryHDU(data=ext.data, header=ext.header).writeto(fixed_file, clobber=True)

        # also save the user-generated position data
        merged = numpy.append(np_all_pixel, np_all_wcs, axis=1)
        numpy.savetxt(datafile, merged)

        ds9.set("frame 3")
        ds9.set("wcs skyformat degrees")
        ds9.set("file %s" % (fixed_file))
        ds9.set("catalog 2mass")
        ds9.set("catalog close")

        ret = "y" #raw_input("happy?")
        if (ret == "r"):
            print "Resetting"
            all_pixel, all_wcs = [], []
            for i, key in enumerate(headers):
                ext.header[key] = backup_hdr[key]
            continue

        happy = (ret.lower() == "y")
    #
    #     next_ext = raw_input("next?")
    #     if (next_ext.lower() == "n"):
    #         break
    #     elif (next_ext.lower() == "wq"):
    #         out_fn = args.filename[:-5]+".fiddlefix.fits"
    #         hdulist.writeto(out_fn, clobber=True)



    #
    # Now we have the rough alignment, get a proper source catalog and
    # re-do the fitting with more accurate positions
    #

    default_param = [
        'X_IMAGE', 'Y_IMAGE',
        'ALPHA_J2000', 'DELTA_J2000',
        'NUMBER',
        'FWHM_IMAGE',
    ]
    with open("tmp.sex.param", "w") as dp:
        dp.write("\n".join(default_param))
    sex_cmd = """sex -CATALOG_NAME tmp.cat -PARAMETERS_NAME tmp.sex.param
        -DETECT_THRESH 10 
        -DETECT_MINAREA 10
        -FILTER N 
        
        %s""" % (fixed_file) # tmp_file
    os.system(" ".join(sex_cmd.split()))

    # load catalog
    catalog = numpy.loadtxt("tmp.cat")

    #
    # Now re-compute the ra/dec of all sources with the rough correction
    # from above - use the fixed file instead
    #


    # match the catalog - ignore all sources with multiple counterparts

    # tree = scipy.spatial.KDTree(catalog[:, 2:4])

    # get ra/dec from reference
    center_pos = wcs.getCentreWCSCoords()
    print center_pos

    coord = astropy.coordinates.SkyCoord(ra=center_pos[0],
                     dec=center_pos[1],
                     unit=(u.degree, u.degree), frame='fk5')
    radius = Quantity(4.0, u.arcmin)
    j = Gaia.cone_search_async(coord, radius)
    r = j.get_results()
    # r.pprint()

    #try:
        #print r[r.keys()[0]]['ra']
    print r['ra']
    print r['dec']
    #except:
    #    pass


    for iteration in range(3):

        # convert the x/y from HST into ra/dec using the latest and greatest WCS
        astro_wcs = astropy.wcs.WCS(
            header=ext.header
        )

        hst_ra_dec = astro_wcs.all_pix2world(catalog[:, 0:2], 1)

        # ra/dec in image
        hst_coords = astropy.coordinates.SkyCoord(
            frame=astropy.coordinates.FK5,
            unit='deg',
            ra=hst_ra_dec[:, 0], dec=hst_ra_dec[:,1],
            # ra=catalog[:, 2], dec=catalog[:, 3],

        )
        print hst_coords


        gaia_pos = astropy.coordinates.SkyCoord(
            frame=astropy.coordinates.FK5,
            unit='deg', ra=r['ra'], dec=r['dec'],
        )
        print gaia_pos

        numpy.savetxt("gaia.cat", numpy.array([r['ra'], r['dec']]).T)


        # result = Vizier.query_region(
        #     coordinates = FK5(ra=Angle(center_pos[0], 'deg'),
        #                       dec=Angle(center_pos[1], 'deg')),
        #     radius=Angle(0.1, "deg"),
        #     catalog='I/324/igsl3', # use GAIA
        #     columns=['RAJ2000','DEJ2000'],
        # )
        # print result
        #
        # print result[result.keys()[0]].pprint()


        matched_cat = astropy.coordinates.match_coordinates_sky(
            matchcoord=hst_coords,
            catalogcoord=gaia_pos,
            nthneighbor=1,
        )
        print matched_cat
        idx, sep2d, dist3d = matched_cat

        good_match = sep2d < Quantity(1.0, u.arcsec)
        print good_match
        print numpy.sum(good_match)

        # TODO: make sure all matches are unique

        print catalog[good_match,0]
        print catalog[good_match,1]

        print r['ra'].shape
        print ((r['ra'])[idx]).shape

        print catalog.shape
        print good_match.shape
        gaia_ra = ((r['ra'])[idx])[good_match]
        print gaia_ra.shape
        gaia_dec = ((r['dec'])[idx])[good_match]

        # gaia_ra = r['ra'][idx][good_match]
        # print gaia_ra
        #
        # matched_gaia = gaia_pos[idx]
        # #print r['ra'][good_match]
        # #print r['dec'][good_match]
        # print matched_gaia
        #
        #
        xy_radec = numpy.array([
             catalog[good_match,0], catalog[good_match,1],
             gaia_ra, gaia_dec,
        ]).T
        print xy_radec
        numpy.savetxt("xy_radec", xy_radec)


        #
        # No recompute the WCS with these new & improved numbers
        #
        fit = scipy.optimize.leastsq(wcs_fit,
                                     p_init,
                                     args=(wcs, headers, xy_radec[:, 0:2], xy_radec[:, 2:4],
                                           cos_dec),
                                     maxfev=1000,
                                     full_output=1)
        # print fit
        final_wcs = fit[0]

        print "WCS correction: %d/%d pixels" % (
            final_wcs[0] - p_init[0],
            final_wcs[1] - p_init[1],
            #     numpy.degrees(improved_wcs[2]-p_init[2]),
        )
        for i,key in enumerate(headers):
            ext.header[key] = final_wcs[i]

    final_file = "tmp_%s.wcsfinal.fits" % (ext.name)
    fits.PrimaryHDU(data=ext.data, header=ext.header).writeto(final_file, clobber=True)
