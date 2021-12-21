use std::path::Path;
use std::fs::File;
use std::io::{BufWriter};
use rand::Rng;
use png;
use vecmath::{self, Vector3, vec3_scale, vec3_add, vec3_sub, vec3_normalized, vec3_dot};

#[derive(Copy,Clone)]
struct Ray {
    a: Vector3<f64>,
    b: Vector3<f64>,
    ht: f64,
    color: Vector3<u8>
}

impl Ray {
    pub fn new(_a: Vector3<f64>, _b: Vector3<f64>, _t: f64) -> Ray {
        Ray { a: _a, b: _b, ht: _t, color: [0,0,0] }
    }

    pub fn origin(self) -> Vector3<f64> {
        self.a
    }
    
    pub fn direction(self) -> Vector3<f64> {
        self.b
    }

    pub fn point_by_t(self, t: f64) -> Vector3<f64> {
        vec3_add(self.a, vec3_scale(self.b, t))
    }
}

struct Camera {
    c: Vector3<f64>,
    tl: Vector3<f64>,
    tr: Vector3<f64>,
    bl: Vector3<f64>,
    br: Vector3<f64>
}

impl Camera {
    pub fn new(_c: Vector3<f64>) -> Camera {
        Camera { 
            c: _c, 
            tl: vec3_add(_c, [-2.0,1.0,-1.0]), 
            tr: vec3_add(_c, [2.0,1.0,-1.0]), 
            bl: vec3_add(_c, [-2.0,-1.0,-1.0]), 
            br: vec3_add(_c, [2.0,-1.0,-1.0]) 
        }
    }
}

#[derive(Copy,Clone)]
struct Sphere {
    c: Vector3<f64>,
    r: f64,
    color: Vector3<u8>
}

impl Sphere {
    pub fn new(_c: Vector3<f64>, _r: f64, _color: Vector3<u8>) -> Sphere {
        Sphere { c: _c, r: _r , color: _color}
    }

    pub fn hit(self, ray: Ray) -> f64 {
        let oc = vec3_sub(ray.origin(), self.c);
        let a = vec3_dot(ray.direction(), ray.direction());
        let b = 2.0 * vec3_dot(oc,ray.direction());
        let c = vec3_dot(oc,oc) - self.r * self.r;
        let discriminant = b*b - 4.0*a*c;

        if discriminant < 0.0 {
            return -1.0;
        }
        else {
            return (-b - discriminant.sqrt()) / (2.0*a);
        }
    }
}

fn main() {
    let path = Path::new(r"test.png");
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let width:usize = 600;
    let height: usize = 300;

    let spheres = vec![Sphere::new([0.0,0.0,-2.0], 1.0, [0,255,0]), Sphere::new([0.0,-800.0,-2.0], 790.0, [0,255,0])];

    let maxt = 1000.0;

    let mut encoder = png::Encoder::new(w, width as u32, height as u32); 
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_trns(vec!(0xFFu8, 0xFFu8, 0xFFu8, 0xFFu8));
    encoder.set_source_gamma(png::ScaledFloat::from_scaled(45455)); 
    encoder.set_source_gamma(png::ScaledFloat::new(1.0 / 2.2));     
    let source_chromaticities = png::SourceChromaticities::new(     
        (0.31270, 0.32900),
        (0.64000, 0.33000),
        (0.30000, 0.60000),
        (0.15000, 0.06000)
    );
    encoder.set_source_chromaticities(source_chromaticities);
    let mut writer = encoder.write_header().unwrap();

    let mut data: [u8; 600 * 300 * 4] = [0; 600 * 300 * 4];
    let mut pixeles = create_rays(Camera::new([0.0,0.0,0.0]), width, height, maxt);
    trace(spheres,&mut pixeles,&mut data, maxt);
    //for i in 0..200 {
        //for j in 0..200 {
            //assign_pixel(i*200+j, ((i as f64/200.0)*255.0) as u8, ((j as f64/200.0)*255.0) as u8, 0, 255, &mut data)
        //}
    //}
    writer.write_image_data(&data).unwrap();
}

fn create_rays(_camera: Camera, _w: usize, _h: usize, t: f64) -> Vec<[Ray;4]>{
    let mut rays: Vec<[Ray;4]> = Vec::new();
    let mut rng = rand::thread_rng();

    let top_left = vec3_sub(_camera.tl,_camera.c);
    let width_offset = vec3_sub(_camera.tr,_camera.tl);
    let height_offset = vec3_sub(_camera.bl,_camera.tl);

    for i in 0.._h {
        for j in 0.._w {
            let n1: f64 = rng.gen();
            let n2: f64 = rng.gen();
            let n3: f64 = rng.gen();
            let n4: f64 = rng.gen();
            let n5: f64 = rng.gen();
            let n6: f64 = rng.gen();
            let n7: f64 = rng.gen();
            let n8: f64 = rng.gen();
            let arr  = [Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n1)/ _h as f64), vec3_scale(width_offset, (j as f64 + n2) / _w as f64))), t),Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n3)/ _h as f64), vec3_scale(width_offset, (j as f64 + n4) / _w as f64))), t),Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n5)/ _h as f64), vec3_scale(width_offset, (j as f64 + n6) / _w as f64))), t),Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n7)/ _h as f64), vec3_scale(width_offset, (j as f64 + n8) / _w as f64))), t)];
            rays.push(arr);
        }
    }

    rays
}


fn trace(spheres: Vec<Sphere>, pixeles: &mut Vec<[Ray;4]>, data: &mut [u8; 600 * 300 * 4], maxt: f64) {
    for i in 0..pixeles.len() {
        let mut pixel_colour: [f64; 3] = [0.0,0.0,0.0];
        for mut j in pixeles[i] {
            for k in 0..spheres.len(){
                let t = spheres[k].hit(j);
                if t > 0.0 && t < j.ht {
                    j.ht = t;
                    let hit_point = j.point_by_t(t);
                    let normal = vec3_scale(vec3_add(vec3_normalized(vec3_sub(hit_point, spheres[k].c)), [1.0,1.0,1.0]),255.0*0.5);
                    j.color = [normal[0] as u8, normal[1] as u8, normal[2] as u8];
                }
            }
            if j.ht == maxt {
                let t = 0.5 * (vec3_normalized(j.direction())[1] + 1.0);
                let c = vec3_add(vec3_scale([255.0,255.0,255.0], 1.0-t),vec3_scale([75.0,150.0,255.0], t));
                j.color = [c[0] as u8, c[1] as u8, c[2] as u8];
            }
            pixel_colour = vec3_add(pixel_colour, [j.color[0] as f64,j.color[1] as f64,j.color[2] as f64]);
        }
        let pc  = vec3_scale(pixel_colour, 0.25);
        assign_pixel(i, pc[0] as u8, pc[1] as u8, pc[2] as u8, 255, data)
    }
}

fn assign_pixel(pixel: usize, r: u8, g: u8, b: u8, a: u8, data: &mut [u8; 600 * 300 * 4]) {
    let i = pixel * 4;
    data[i] = r;
    data[i+1] = g;
    data[i+2] = b;
    data[i+3] = a;
}
