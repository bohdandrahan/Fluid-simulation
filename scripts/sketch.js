let N = 256
let dt = 0.1;
let diffusion = 0;
let viscosity = 0;
let fluid;
let iter = 1;

function IX(x, y) {
	x = constrain(x, 0, N - 1)
	y = constrain(y, 0, N - 1)
	return x + y * N
}

function setup() {
	createCanvas(N, N)
	fluid = new Fluid(dt, diffusion, viscosity);

}

function mouseDragged() {
	fluid.addDye(mouseX, mouseY, 100)
	let amtX = mouseX - pmouseX;
	let amtY = mouseY - pmouseY;
	fluid.addVelocity(mouseX, mouseY, amtX, amtY)

}

function draw() {
	background(100)
	fluid.step()
	fluid.renderD()
	// put drawing code here
}