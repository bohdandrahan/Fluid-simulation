let N = 100
let dt = 0.1;
let diffusion = 0;
let viscosity = 0;
let fluid;
let iter = 1;
let scale = 10

function IX(x, y) {
	x = constrain(x, 0, N - 1)
	y = constrain(y, 0, N - 1)
	return x + y * N
}

function setup() {
	createCanvas(N * scale, N * scale)
	fluid = new Fluid(dt, diffusion, viscosity);

}

function mouseDragged() {
	fluid.addDye(floor(mouseX / scale), floor(mouseY / scale), 100)
	let amtX = random(0, 1);
	let amtY = random(0, 1);
	fluid.addVelocity(floor(mouseX / scale), floor(mouseY / scale), amtX, amtY)

}

function draw() {
	background(100)
	fluid.step()
	fluid.renderD()
	// put drawing code here
}