class Fluid {
	constructor(dt, diffusion, viscosity) {
		this.size = N;
		this.dt = dt;
		this.diffusion = diffusion;
		this.viscosity = viscosity;


		this.s = new Array(N * N);
		for (let i = 0; i < N * N; i++) {
			this.s[i] = 0
		}
		this.density = new Array(N * N);
		for (let i = 0; i < N * N; i++) {
			this.density[i] = 0
		}
		this.Vx = new Array(N * N);
		for (let i = 0; i < N * N; i++) {
			this.Vx[i] = 0
		}
		this.Vy = new Array(N * N);
		for (let i = 0; i < N * N; i++) {
			this.Vy[i] = 0
		}

		this.Vx0 = new Array(N * N);
		for (let i = 0; i < N * N; i++) {
			this.Vx0[i] = 0
		}
		this.Vy0 = new Array(N * N);
		for (let i = 0; i < N * N; i++) {
			this.Vy0[i] = 0
		}
	};
	addDye(x, y, amount) {
		let index = IX(x, y);
		this.density[index] += amount;
	}
	step() {

		let visc = this.viscosity;
		let diff = this.diffusion;
		let dt = this.dt;
		let Vx = this.Vx;
		let Vy = this.Vy;
		let Vx0 = this.Vx0;
		let Vy0 = this.Vy0;
		let s = this.s;
		let density = this.density;


		diffuse(1, Vx0, Vx, visc, dt);
		diffuse(2, Vy0, Vy, visc, dt);
		project(Vx0, Vy0, Vx, Vy);

		advect(1, Vx, Vx0, Vx0, Vy0, dt);
		advect(2, Vy, Vy0, Vx0, Vy0, dt);

		project(Vx, Vy, Vx0, Vy0);

		diffuse(0, s, density, diff, dt);
		advect(0, density, s, Vx, Vy, dt);

	}
	renderD() {
		for (let i = 0; i < N; i++) {
			for (let j = 0; j < N; j++) {
				let d = this.density[IX(i, j)]
				let x = i * scale;
				let y = j * scale;
				fill(d)
				stroke(100, 0, 0)
				rect(x, y, scale, scale)
			}
		}
	}

	addVelocity(x, y, dVx, dVy) {
		let index = IX(x, y);
		this.Vx[index] += dVx;
		this.Vy[index] += dVy;
	}
}

function diffuse(b, x, x0, diffusion, dt) {
	let a = dt * diffusion * (N - 2) * (N - 2)
	lin_solve(b, x, x0, a, 1 + 6 * a);
}

function lin_solve(b, x, x0, a, c) {
	let cRecip = 1.0 / c;
	for (let j = 1; j < N - 1; j++) {
		for (let i = 1; i < N - 1; i++) {
			x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip
		}
	}

}

function project(velocX, velocY, p, div) {
	for (let j = 1; j < N - 1; j++) {
		for (let i = 1; i < N - 1; i++) {
			div[IX(i, j)] = -0.5 * (
				velocX[IX(i + 1, j)] -
				velocX[IX(i - 1, j)] +
				velocY[IX(i, j + 1)] -
				velocY[IX(i, j - 1)]
			) / N;
			p[IX(i, j)] = 0;
		}

	}
	set_bnd(0, div);
	set_bnd(0, p);
	lin_solve(0, p, div, 1, 6);

	for (let j = 1; j < N - 1; j++) {
		for (let i = 1; i < N - 1; i++) {
			velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
			velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N - p[IX(i, j)] * N
		}
	}



	set_bnd(1, velocX);
	set_bnd(2, velocY);
}


function advect(b, d, d0, velocX, velocY, dt) {
	let i0, i1, j0, j1;

	let dtx = dt * (N - 2);
	let dty = dt * (N - 2);

	let s0, s1, t0, t1;
	let tmp1, tmp2, tmp3, x, y;

	let Nfloat = N;
	let ifloat, jfloat, kfloat;
	let i, j;

	for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
		for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
			tmp1 = dtx * velocX[IX(i, j)];
			tmp2 = dty * velocY[IX(i, j)];
			x = ifloat - tmp1;
			y = jfloat - tmp2;

			if (x < 0.5) {
				x = 0.5
			};
			if (x > Nfloat + 0.5) {
				x = Nfloat + 0.5
			};
			i0 = floor(x);
			i1 = i0 + 1.0;
			if (y < 0.5) {
				y = 0.5
			};
			if (y > Nfloat + 0.5) {
				y = Nfloat + 0.5
			};
			j0 = floor(y);
			j1 = j0 + 1.0;

			s1 = x - i0;
			s0 = 1.0 - s1;
			t1 = y - j0;
			t0 = 1.0 - t1;

			let i0i = i0;
			let i1i = i1;
			let j0i = j0;
			let j1i = j1;

			d[IX(i, j)] = s0 * t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)] + s1 * t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)];
		}
	}
	set_bnd(b, d);
}

function set_bnd(b, x, N) {

	for (let i = 1; i < N - 1; i++) {
		x[IX(i, 0)] = b === 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N - 1)] = b === 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
	}

	for (let j = 1; j < N - 1; j++) {
		x[IX(0, j)] = b === 1 ? -x[IX(1, j)] : x[IX(1, j)];
		x[IX(N - 1, j)] = b === 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
	}
	x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N - 1)] = 0.5 * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
	x[IX(N - 1, 0)] = 0.5 * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
	x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);

}