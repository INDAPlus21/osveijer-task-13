use std::path::Path;
use std::fs::File;
use std::io::{BufWriter};
use rand::Rng;
use png;
use vecmath::{self, Vector3, vec3_scale, vec3_add, vec3_sub, vec3_normalized, vec3_dot, vec3_len};

#[derive(Copy,Clone)]
struct Ray {
    a: Vector3<f64>,
    b: Vector3<f64>,
}

impl Ray {
    pub fn new(_a: Vector3<f64>, _b: Vector3<f64>) -> Ray {
        Ray { a: _a, b: _b}
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

    pub fn hit(self, ray: Ray, mint: f64, maxt: f64, rec: &mut HitRecord) -> bool {
        let oc = vec3_sub(ray.origin(), self.c);
        let a = vec3_dot(ray.direction(), ray.direction());
        let b = 2.0 * vec3_dot(oc,ray.direction());
        let c = vec3_dot(oc,oc) - self.r * self.r;
        let discriminant = b*b - 4.0*a*c;

        if discriminant > 0.0 {
            let t1 = (-b - discriminant.sqrt()) / (2.0*a);
            if t1 < maxt && t1 > mint {
                rec.t = t1;
                rec.p = ray.point_by_t(t1);
                rec.normal = vec3_normalized(vec3_sub(rec.p, self.c));
                return true;
            }
            let t2 = (-b + discriminant.sqrt()) / (2.0*a);
            if t2 < maxt && t2 > mint {
                rec.t = t2;
                rec.p = ray.point_by_t(t2);
                rec.normal = vec3_normalized(vec3_sub(rec.p, self.c));
                return true;
            }
        }
        false
    }
}

struct HitRecord {
    t: f64,
    p: Vector3<f64>,
    normal: Vector3<f64>
}

impl HitRecord {
    pub fn new() -> HitRecord {
        HitRecord {
            t: 0.0,
            p: [0.0,0.0,0.0],
            normal: [0.0,0.0,0.0]
        }
    }
}

fn main() {
    let path = Path::new(r"test.png");
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let width:usize = 600;
    let height: usize = 300;

    let spheres = vec![Sphere::new([0.0,0.0,-2.0], 1.0, [0,255,0]), Sphere::new([0.0,-800.0,-2.0], 799.0, [0,255,0])];

    let maxt = 1000.0;
    let mint = 0.0001;

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
    trace_all(spheres,&mut pixeles,&mut data, mint, maxt);
    
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
            let arr  = [Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n1)/ _h as f64), vec3_scale(width_offset, (j as f64 + n2) / _w as f64)))),Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n3)/ _h as f64), vec3_scale(width_offset, (j as f64 + n4) / _w as f64)))),Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n5)/ _h as f64), vec3_scale(width_offset, (j as f64 + n6) / _w as f64)))),Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, (i as f64 + n7)/ _h as f64), vec3_scale(width_offset, (j as f64 + n8) / _w as f64))))];
            rays.push(arr);
        }
    }

    rays
}


fn trace_all(spheres: Vec<Sphere>, pixeles: &mut Vec<[Ray;4]>, data: &mut [u8; 600 * 300 * 4], mint: f64, maxt: f64) {
    for i in 0..pixeles.len() {
        let mut pixel_colour: [f64; 3] = [0.0,0.0,0.0];
        for j in pixeles[i] {
            pixel_colour = vec3_add(pixel_colour, trace(j, &spheres, mint, maxt));
        }
        let mut pc  = vec3_scale(pixel_colour, 0.25);
        pc = [pc[0].sqrt(), pc[1].sqrt(),pc[2].sqrt()];
        pc = vec3_scale(pc, 255.0);
        assign_pixel(i, pc[0] as u8, pc[1] as u8, pc[2] as u8, 255, data)
    }
}

fn trace(ray: Ray, spheres: &Vec<Sphere>, mint: f64, maxt: f64) -> Vector3<f64> {
    let mut rec: HitRecord = HitRecord::new();
    if hits(ray, spheres, mint, maxt, &mut rec) {
        let new_ray_dir = vec3_sub(vec3_add(vec3_add(rec.p, rec.normal), get_random_in_sphere()), rec.p);
        return vec3_scale(trace(Ray::new(rec.p,new_ray_dir),spheres,mint,maxt), 0.5);
    }
    else {
        let t = 0.5 * (vec3_normalized(ray.direction())[1] + 1.0);
        return vec3_add(vec3_scale([1.0,1.0,1.0], 1.0-t),vec3_scale([0.3,0.5,1.0], t));
    }
}

fn hits(ray: Ray, spheres: &Vec<Sphere>, mint: f64, maxt: f64, rec: &mut HitRecord) -> bool {
    let mut temp_rec: HitRecord = HitRecord::new();
    let mut hit_anything = false;
    let mut min_dist = maxt;
    for k in 0..spheres.len(){
        if spheres[k].hit(ray,mint,min_dist,&mut temp_rec) {
            hit_anything = true;
            min_dist = temp_rec.t;
            rec.t = temp_rec.t;
            rec.p = temp_rec.p;
            rec.normal = temp_rec.normal;
        }
    }
    hit_anything
}

fn get_random_in_sphere() -> Vector3<f64> {
    let mut rng = rand::thread_rng();
    let mut vec: Vector3<f64> = [2.0,2.0,2.0];
    while vec3_len(vec) > 1.0 {
        let x: f64 = rng.gen();
        let y: f64 = rng.gen();
        let z: f64 = rng.gen();
        vec = vec3_sub(vec3_scale([x,y,z], 2.0),[1.0,1.0,1.0]);
    }
    vec
}

fn assign_pixel(pixel: usize, r: u8, g: u8, b: u8, a: u8, data: &mut [u8; 600 * 300 * 4]) {
    let i = pixel * 4;
    data[i] = r;
    data[i+1] = g;
    data[i+2] = b;
    data[i+3] = a;
}
