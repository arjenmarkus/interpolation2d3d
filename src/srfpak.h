/* srfpak.h --
       Header file for the Fortran-C interface to the SRFPAK package
*/

/* Opaque type for storing the network data */
typedef struct {
    void *network_f;
} network_2d;

/* createMesh --
       Create the mesh (and store the result in the network struct)

   Arguments:
       network           The network struct to be used
       n                 Number of points (network nodes)
       x                 X-coordinates of the points
       y                 Y-coordinates of the points
*/
void createMesh( network_2d *network, int n, float *x, float *y );

/* deleteMesh --
       Delete the contents of the network struct

   Arguments:
       network           The network struct to be freed up
*/
void deleteMesh( network_2d *network );

/* setMeshValues --
       Set the values at the network nodes (points)

   Arguments:
       network           The network struct to be used
       z                 Values at the points
*/
void setMeshValues( network_2d *network, float *z );

/* interpolatePoint --
       Interpolate the data for a single location (px,py)

   Arguments:
       network           The network struct to be used
       px                X-coordinate of the point
       py                Y-coordinate of the point
*/
float interpolatePoint( network_2d *network, float x, float y );

/* meshVolume --
       Volume as occupied by the points and the values in the network

   Arguments:
       network           The network struct to be used
*/
float meshVolume( network_2d *network );

