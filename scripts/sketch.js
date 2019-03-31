let N = 256
let dt = 0.1;
let diffusion = 0;
let viscosity = 0;
let fluid;
let iter = 10;

function IX(x, y) {
	return x + y * N
}

function setup() {
	createCanvas(N, N)
	fluid = new Fluid(dt, diffusion, viscosity);
}

function draw() {
	background(100)
	// put drawing code here
}