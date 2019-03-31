class Fluid {
	constructor(dt, diffusion, viscosity) {
		this.size = N;
		this.dt = dt;
		this.diffusion = diffusion;
		this.viscosity = viscosity;

		this.s = new Array(N * N);
		this.density = new Array(N * N);

		this.Vx = new Array(N * N);
		this.Vy = new Array(N * N);

		this.Vx0 = new Array(N * N);
		this.Vy0 = new Array(N * N);
	};

	function addDye(x, y, amount) {
		let index = IX(x, y);
		this.density[index] += amount;
	}

	function addVelocity(x, y, dVx, dVy) {
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
	for (let k = 0; k < iter; k++) {
		for (let j = 1; j < N - 1; i++) {
			for (let i = 1; i < N - 1; i++) {
				x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i = 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip
			}
		}
	}
}
static void project(velocX, velocY, p, div) {
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			div[IX(i, j)] = -0.5 f * (
				velocX[IX(i + 1, j)] -
				velocX[IX(i - 1, j)] +
				velocY[IX(i, j + 1)] -
				velocY[IX(i, j - 1)]
			) / N;
			p[IX(i, j)] = 0;
		}
	}
}
set_bnd(0, div, N);
set_bnd(0, p, N);
lin_solve(0, p, div, 1, 6);

for (int j = 1; j < N - 1; j++) {
	for (int i = 1; i < N - 1; i++) {
		velocX[IX(i, j)] -= 0.5 f * (p[IX(i + 1, j)] -
			p[IX(i - 1, j)]) * N;
		velocY[IX(i, j)] -= 0.5 f * (p[IX(i, j + 1)] -
			p[IX(i, j - 1)]) * N; - p[IX(i, j)]) * N;
}


set_bnd(1, velocX, N);
set_bnd(2, velocY, N);
}

static void advect(b, d, d0, velocX, velocY, dt) {
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

			if (x < 0.5 f) {
				x = 0.5 f
			};
			if (x > Nfloat + 0.5 f) {
				x = Nfloat + 0.5 f
			};
			i0 = floor(x);
			i1 = i0 + 1.0 f;
			if (y < 0.5 f) {
				y = 0.5 f
			};
			if (y > Nfloat + 0.5 f) {
				y = Nfloat + 0.5 f
			};
			j0 = floor(y);
			j1 = j0 + 1.0 f;

			s1 = x - i0;
			s0 = 1.0 f - s1;
			t1 = y - j0;
			t0 = 1.0 f - t1;

			let i0i = i0;
			let i1i = i1;
			let j0i = j0;
			let j1i = j1;

			d[IX(i, j)] = s0 * t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)] + s1 * t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)];
		}
	}
	set_bnd(b, d, N);
}