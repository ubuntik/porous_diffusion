/*
 *
 */

void get_K(matrix &K)
{
	K(0, 0) = 700.72; K(1, 0) = -90.791; K(2, 0) = -51.8508; K(3, 0) = -20.9356;
	K(4, 0) = -16.6301; K(5, 0) = -36.7135; K(6, 0) = -12.6038; K(7, 0) = -11.0965;

	K(0, 1) = -90.4545; K(1, 1) = 410.508; K(2, 1) = 10.7879; K(3, 1) = -6.02379;
	K(4, 1) = 7.86314; K(5, 1) = 4.83955; K(6, 1) = 13.3205; K(7, 1) = 6.58391;

	K(0, 2) = -51.7445; K(1, 2) = 10.9115; K(2, 2) = 162.8; K(3, 2) = 24.564;
	K(4, 2) = 12.9702; K(5, 2) = 16.7895; K(6, 2) = 16.4923; K(7, 2) = 20.6626;

	K(0, 3) = -20.952; K(1, 3) = -5.76781; K(2, 3) = 24.6143; K(3, 3) = 110.027;
	K(4, 3) = 7.23522; K(5, 3) = 14.8657; K(6, 3) = 10.6427; K(7, 3) = 18.3762;

	K(0, 4) = -16.8055; K(1, 4) = 7.7814; K(2, 4) = 13.0678; K(3, 4) = 7.29436;
	K(4, 4) = 85.9633; K(5, 4) = 10.5721; K(6, 4) = 10.8133; K(7, 4) = 13.1033;

	K(0, 5) = -36.7031; K(1, 5) = 4.94953; K(2, 5) = 16.7641; K(3, 5) = 14.8475;
	K(4, 5) = 10.5869; K(5, 5) = 91.9372; K(6, 5) = 13.9441; K(7, 5) = 18.5697;

	K(0, 6) = -12.6047; K(1, 6) = 13.253; K(2, 6) = 16.6121; K(3, 6) = 10.6599;
	K(4, 6) = 10.821; K(5, 6) = 13.9652; K(6, 6) = 85.7859; K(7, 6) = 23.7649;

	K(0, 7) = -11.1; K(1, 7) = 6.64; K(2, 7) = 20.7; K(3, 7) = 18.4;
	K(4, 7) = 13.1; K(5, 7) = 18.6; K(6, 7) = 23.8; K(7, 7) = 68.8;
}

static void get_D(matrix &D)
{
	D(0, 0) = 0.00579; D(1, 0) = -0.00135; D(2, 0) = -0.00060; D(3, 0) = -0.00051;
	D(4, 0) = -0.00054; D(5, 0) = -0.00072; D(6, 0) = -0.00097; D(7, 0) = -0.00111;

	D(0, 1) = -0.00135; D(1, 1) = 0.00462; D(2, 1) = -0.00040; D(3, 1) = -0.00038;
	D(4, 1) = -0.00041; D(5, 1) = -0.00053; D(6, 1) = -0.00070; D(7, 1) = -0.00083;

	D(0, 2) = -0.00060; D(1, 2) = -0.00040; D(2, 2) = 0.00242; D(3, 2) = -0.00018;
	D(4, 2) = -0.00019; D(5, 2) = -0.00028; D(6, 2) = -0.00035; D(7, 2) = -0.00043;

	D(0, 3) = -0.00051; D(1, 3) = -0.00038; D(2, 3) = -0.00018; D(3, 3) = 0.00219;
	D(4, 3) = -0.00017; D(5, 3) = -0.00025; D(6, 3) = -0.00031; D(7, 3) = -0.00039;

	D(0, 4) = -0.00054; D(1, 4) = -0.00041; D(2, 4) = -0.00019; D(3, 4) = -0.00017;
	D(4, 4) = 0.00228; D(5, 4) = -0.00024; D(6, 4) = -0.00033; D(7, 4) = -0.00041;

	D(0, 5) = -0.00072; D(1, 5) = -0.00053; D(2, 5) = -0.00028; D(3, 5) = -0.00025;
	D(4, 5) = -0.00024; D(5, 5) = 0.00304; D(6, 5) = -0.00045; D(7, 5) = -0.00056;

	D(0, 6) = -0.00097; D(1, 6) = -0.00070; D(2, 6) = -0.00035; D(3, 6) = -0.00031;
	D(4, 6) = -0.00033; D(5, 6) = -0.00045; D(6, 6) = 0.00383; D(7, 6) = -0.00071;

	D(0, 7) = -0.00112; D(1, 7) = -0.00083; D(2, 7) = -0.00043; D(3, 7) = -0.00039;
	D(4, 7) = -0.00041; D(5, 7) = -0.00056; D(6, 7) = -0.00071; D(7, 7) = 0.00445;
}

