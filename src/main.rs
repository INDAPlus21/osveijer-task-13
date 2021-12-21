use std::path::Path;
use std::fs::File;
use std::io::{BufWriter};
use png;
use vecmath::{self, Vector3, vec3_scale, vec3_add, vec3_sub, vec3_normalized, vec3_dot};

#[derive(Copy,Clone)]
struct Ray {
    a: Vector3<f32>,
    b: Vector3<f32>
}

impl Ray {
    pub fn new(_a: Vector3<f32>, _b: Vector3<f32>) -> Ray {
        Ray { a: _a, b: _b }
    }

    pub fn origin(self) -> Vector3<f32> {
        self.a
    }
    
    pub fn direction(self) -> Vector3<f32> {
        self.b
    }

    pub fn point_by_t(self, t: f32) -> Vector3<f32> {
        vec3_add(self.a, vec3_scale(self.b, t))
    }
}

struct Camera {
    c: Vector3<f32>,
    tl: Vector3<f32>,
    tr: Vector3<f32>,
    bl: Vector3<f32>,
    br: Vector3<f32>
}

impl Camera {
    pub fn new(_c: Vector3<f32>) -> Camera {
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
    c: Vector3<f32>,
    r: f32,
    color: Vector3<u8>
}

impl Sphere {
    pub fn new(_c: Vector3<f32>, _r: f32, _color: Vector3<u8>) -> Sphere {
        Sphere { c: _c, r: _r , color: _color}
    }

    pub fn hit(self, ray: Ray) -> bool {
        let oc = vec3_sub(ray.origin(), self.c);
        let a = vec3_dot(ray.direction(), ray.direction());
        let b = 2.0 * vec3_dot(oc,ray.direction());
        let c = vec3_dot(oc,oc) - self.r * self.r;
        let discriminant = b*b - 4.0*a*c;

        return discriminant > 0.0;
    }
}

fn main() {
    let path = Path::new(r"test.png");
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let width:usize = 200;
    let height: usize = 100;

    let spheres = vec![Sphere::new([0.0,0.0,-2.0], 1.0, [0,255,0])];

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

    let mut data: [u8; 200 * 100 * 4] = [0; 200 * 100 * 4];
    let rays = create_rays(Camera::new([0.0,0.0,0.0]), width, height);
    set_background(&rays, &mut data);
    check_spheres(spheres,&rays,&mut data);
    //for i in 0..200 {
        //for j in 0..200 {
            //assign_pixel(i*200+j, ((i as f32/200.0)*255.0) as u8, ((j as f32/200.0)*255.0) as u8, 0, 255, &mut data)
        //}
    //}
    writer.write_image_data(&data).unwrap();
}

fn create_rays(_camera: Camera, _w: usize, _h: usize) -> Vec<Ray>{
    let mut rays: Vec<Ray> = Vec::new();

    let top_left = vec3_sub(_camera.tl,_camera.c);
    let width_offset = vec3_sub(_camera.tr,_camera.tl);
    let height_offset = vec3_sub(_camera.bl,_camera.tl);

    for i in 0.._h {
        for j in 0.._w {
            rays.push(Ray::new(_camera.c, vec3_add(top_left,vec3_add(vec3_scale(height_offset, i as f32 / _h as f32), vec3_scale(width_offset, j as f32 / _w as f32)))));
        }
    }

    rays
}

fn set_background(rays: &Vec<Ray>, data: &mut [u8; 200 * 100 * 4]) {
    for i in 0..(rays.len() - 1) {
        let t = 0.5 * (vec3_normalized(rays[i].direction())[1] + 1.0);
        let c = vec3_add(vec3_scale([255.0,255.0,255.0], 1.0-t),vec3_scale([75.0,150.0,255.0], t));
        assign_pixel(i, c[0] as u8, c[1] as u8, c[2] as u8, 255, data);
    }
}

fn check_spheres(spheres: Vec<Sphere>, rays: &Vec<Ray>, data: &mut [u8; 200 * 100 * 4]) {
    for i in spheres {
        for j in 0..(rays.len() - 1) {
            if i.hit(rays[j]) {
                assign_pixel(j, i.color[0], i.color[1], i.color[2], 255, data);
            }
        }
    }
}

fn assign_pixel(pixel: usize, r: u8, g: u8, b: u8, a: u8, data: &mut [u8; 200 * 100 * 4]) {
    let i = pixel * 4;
    data[i] = r;
    data[i+1] = g;
    data[i+2] = b;
    data[i+3] = a;
}
