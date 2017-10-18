package main

import (
	"fmt"
	"image"
	"log"
	"math"
	"math/cmplx"
	"os"

	"github.com/gonum/matrix/mat64"
	"github.com/mjibson/go-dsp/fft"

	"image/color"
	"image/png"
)

// "fmt"
// "image"
// "log"
// "os"
//
// "image/color"
// _ "image/gif"
// _ "image/jpeg"

func main() {

	// c
	f_c, err := os.Open("c.png")
	if err != nil {
		log.Fatal(err)
	}
	defer f_c.Close()
	im_c, _, err := image.Decode(f_c)
	if err != nil {
		log.Fatal(err)
	}
	bounds := im_c.Bounds()
	width, height := bounds.Max.X-bounds.Min.X, bounds.Max.Y-bounds.Min.Y
	h_c := make([]float64, width*height)
	hi_c := make([][]complex128, width)
	for c := 0; c < width; c++ {
		hi_c[c] = make([]complex128, height)
		for r := 0; r < height; r++ {
			h_c[c*height+r] = float64(color.GrayModel.Convert(im_c.At(c, r)).(color.Gray).Y)
			hi_c[c][r] = complex(h_c[c*height+r], 0)
		}
	}
	fmt.Printf("read %s, size [%dx%d]\n", f_c.Name(), width, height)

	// n
	h_n := make([]float64, width*height)
	a := 2.0
	a2 := a * a
	s := 0.0
	dw := 5
	dh := 5
	for c := -dw; c < dw; c++ {
		for r := -dh; r < dh; r++ {
			// apparatus function is around 0,0 point
			nc := (c + width) % width
			nr := (r + height) % height

			x2 := float64(nc * nc)
			y2 := float64(nr * nr)
			h_n[nc*height+nr] = math.Exp((-x2-y2)/a2) / (a2 * math.Pi)

			s += h_n[nc*height+nr]
		}
	}
	hi_n := make([][]complex128, width)
	for c := 0; c < width; c++ {
		hi_n[c] = make([]complex128, height)
	}
	for c := -dw; c < dw; c++ {
		for r := -dh; r < dh; r++ {
			// apparatus function is around 0,0 point
			nc := (c + width) % width
			nr := (r + height) % height

			h_n[nc*height+nr] /= s
			hi_n[nc][nr] = complex(h_n[nc*height+nr], 0)
		}
	}

	// FFT
	fft_c := fft.FFT2(hi_c)
	fft_n := fft.FFT2(hi_n)

	abs_fft_c, h_max := abs(fft_c, width, height)
	fmt.Printf("h_max = %f\n", h_max)

	// save fft_c
	im_fft_c := image.NewGray(image.Rect(0, 0, width, height))
	for c := 0; c < width; c++ {
		for r := 0; r < height; r++ {
			h := math.Log1p(abs_fft_c[c*height+r]) * 255 / math.Log1p(h_max)
			im_fft_c.Set(c, r, color.Gray{uint8(h)})
		}
	}
	f_fft_c, err := os.Create("fft_c.png")
	defer f_fft_c.Close()
	png.Encode(f_fft_c, im_fft_c)
	fmt.Printf("saved %s\n", f_fft_c.Name())

	// for r := 2; r < 200; r++ {
	// 	_, _, _, deltaI := approx(abs_fft_c, width, height, r, -1)
	// 	xaE, xbE, xcE, deltaE := approx(abs_fft_c, width, height, r, 1)
	// 	fmt.Printf("%d: %f	%f	%f	%f	%f	%f\n", r, xaE, xbE, xcE, deltaE, deltaI, deltaE+deltaI)
	// }

	//!xaE, xbE, xcE, _ := approx(fft_c, width, height, 50, 1)

	// un_c
	fft_un_c := make([][]complex128, width)
	k_pos := 0
	k_neg := 0
	for c := 0; c < width; c++ {
		fft_un_c[c] = make([]complex128, height)
		for r := 0; r < height; r++ {
			k := 0.0

			//!x := float64(c)
			//!y := float64(r)
			//!if x > float64(width/2) {
			//!	x = float64(width) - x
			//!}
			//!if y > float64(height/2) {
			//!	y = float64(height) - y
			//!}
			//!N2 := (xaE*x + xbE*y + xcE) * (xaE*x + xbE*y + xcE)
			//!C2 := cmplx.Abs(fft_c[c][r]) * cmplx.Abs(fft_c[c][r])
			//!
			//!k = 1.0 - N2/C2
			//!if k > 0 {
			//!	k_pos += 1
			//!} else {
			//!	k = 0
			//!	k_neg += 1
			//!}

			// just ckeck
			//if x*x+y*y > 100 {
			cc := c //(c + width/2) % width
			rr := r //(r + height/2) % height
			if cc*cc+rr*rr > 100 {
				k = 0
			} else {
				k = 1
			}

			Z := fft_c[c][r] / fft_n[c][r]
			fft_un_c[c][r] = complex(real(Z)*k, imag(Z))
		}
	}
	fmt.Printf("k_pos %d, k_neg %d\n", k_pos, k_neg)
	un_c := fft.IFFT2(fft_un_c)

	// save un_c
	im_un_c := image.NewGray(image.Rect(0, 0, width, height))
	for c := 0; c < width; c++ {
		for r := 0; r < height; r++ {
			height := cmplx.Abs(un_c[c][r])
			if height > 255 {
				height = 255
			}
			im_un_c.Set(c, r, color.Gray{uint8(height)})
		}
	}
	f_un_c, err := os.Create("un_c.png")
	defer f_un_c.Close()
	png.Encode(f_un_c, im_un_c)
	fmt.Printf("saved %s\n", f_un_c.Name())
}

// side = 1 => exterior
// side = -1 => interior
//
// img in [0..width/2, 0..height/2]
func approx(img [][]complex128, width int, height int, dr int, side int) (xa, xb, xc, delta float64) {
	xa = 0
	xb = 0
	xc = 0
	delta = 0
	// (ax^2+by^2+cxy+d-f(x,y))*x^2 = 0
	//                          y^2 = 0
	//                          xy  = 0
	//                          1   = 0
	//
	// x4   x2y2 x3y  x2  | fx2
	// x2y2 y4   xy3  y2  | fy2
	// x3y  xy3  x2y2 xy  | fxy
	// x2   y2   xy   1   | f

	// instead use a plate:

	// (ax + by + c - f(x,y))^2 ~> min
	//
	// (ax + by + c - f(x,y))x
	// (ax + by + c - f(x,y))y
	// (ax + by + c - f(x,y))
	//
	// x2 xy x1 | fx
	// xy y2 y1 | fy
	// x1 y1  1 | f

	n := 0
	for c := 0; c < width/2; c++ {
		for r := 0; r < height/2; r++ {
			if side*(c*c+r*r) > side*(dr*dr) {
				n = n + 1
			}
		}
	}

	var (
		x2, y2, xy, x1, y1, I float64
		F, Fx, Fy             float64
	)
	dd := 1.0 / float64(n)
	for c := 0; c < width/2; c++ {
		for r := 0; r < height/2; r++ {
			if side*(c*c+r*r) > side*(dr*dr) {
				f_im := cmplx.Abs(img[c][r])

				x := c
				y := r

				x2 += float64(x*x) * dd
				y2 += float64(y*y) * dd
				xy += float64(x*y) * dd
				x1 += float64(x) * dd
				y1 += float64(y) * dd
				I += dd

				F += f_im * dd
				Fx += float64(x) * f_im * dd
				Fy += float64(y) * f_im * dd
			}
		}
	}

	a := mat64.NewDense(3, 3, []float64{
		x2, xy, x1,
		xy, y2, y1,
		x1, y1, I,
	})
	b := mat64.NewVector(3, []float64{Fx, Fy, F})

	var x mat64.Vector
	if err := x.SolveVec(a, b); err != nil {
		//fmt.Println("Matrix is near singular: ", err)
	}
	//fmt.Printf("x = %0.4v\n", mat64.Formatted(&x, mat64.Prefix("    ")))

	xa = x.At(0, 0)
	xb = x.At(1, 0)
	xc = x.At(2, 0)
	for c := 0; c < width/2; c++ {
		for r := 0; r < height/2; r++ {
			if side*(c*c+r*r) > side*(dr*dr) {
				f_im := cmplx.Abs(img[c][r])

				x := float64(c)
				y := float64(r)

				delta += math.Sqrt((xa*x+xb*y+xc-f_im)*(xa*x+xb*y+xc-f_im)) * dd
			}
		}
	}

	return
}

// img in [0..width/2, 0..height/2]
func abs(img [][]complex128, width, height int) ([]float64, float64) {
	h := make([]float64, width*height)
	h_max := 0.0
	for c := 0; c < width; c++ {
		for r := 0; r < height; r++ {
			cc := (c + width/2) % width
			rr := (r + height/2) % height

			h[cc*height+rr] = cmplx.Abs(img[c][r])
			if h[cc*height+rr] > h_max {
				h_max = h[cc*height+rr]
			}
		}
	}
	return h, h_max
}

// func dump(img []float64, width, height int, nm string) {
// 	f, err := os.Create(nm + ".pgm")
// 	if err != nil {
// 		log.Fatal(err)
// 	}
// 	defer f.Close()
// 	png.Encode(f_un_c, im_un_c)
// 	fmt.Printf("saved %s\n", f_un_c.Name())
// }
