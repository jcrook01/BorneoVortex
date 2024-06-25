import warnings
import iris

def set_latlon_bounds(cube):
    for axis in 'XY':
        coord = cube.coord(axis=axis)
        if len(coord.points)==1:
            return
        # If no bounds, set them
        if not coord.has_bounds():
            coord.guess_bounds()
        # If bounds are not contiguous, get rid of them and make them
        # again
        elif not coord.is_contiguous():
            coord.bounds = None
            coord.guess_bounds()
        # If bounds are still not contiguous, adjust them
        if not coord.is_contiguous():
            b = coord._bounds.copy()
            for (ii, (b1, b2)) in enumerate(zip(b[:-1], b[1:])):
                if b1[1] != b2[0]:
                    b[ii][1] = b[ii+1][0]
            coord._bounds = b

if __name__=='__main__':
    main()
