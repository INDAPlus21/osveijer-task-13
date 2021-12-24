use std::path::Path;
use std::fs::File;
use std::io::{BufWriter};
use rand::Rng;
use png;
use vecmath::{self, Vector3, vec3_scale, vec3_add, vec3_sub, vec3_normalized, vec3_dot, vec3_len, vec3_mul, vec3_cross};
use rayon::prelude::*;

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
}

impl Camera {
    pub fn new(_c: Vector3<f64>) -> Camera {
        Camera { 
            c: _c, 
            tl: vec3_add(_c, [-2.0,1.0,-1.0]), 
            tr: vec3_add(_c, [2.0,1.0,-1.0]), 
            bl: vec3_add(_c, [-2.0,-1.0,-1.0]), 
        }
    }
}

#[derive(Copy,Clone)]
struct Sphere {
    c: Vector3<f64>,
    r: f64,
    mat: Material
}

impl Sphere {
    pub fn new(_c: Vector3<f64>, _r: f64, _mat: Material) -> Sphere {
        Sphere { c: _c, r: _r , mat: _mat}
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
                rec.mat = self.mat;
                return true;
            }
            let t2 = (-b + discriminant.sqrt()) / (2.0*a);
            if t2 < maxt && t2 > mint {
                rec.t = t2;
                rec.p = ray.point_by_t(t2);
                rec.normal = vec3_normalized(vec3_sub(rec.p, self.c));
                rec.mat = self.mat;
                return true;
            }
        }
        false
    }
}

#[derive(Copy,Clone)]
struct Plane {
    normal: Vector3<f64>,
    holder: f64,
    center: Vector3<f64>,
    mat: Material,
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
    min_z: f64,
    max_z: f64
}

impl Plane {
    pub fn new(_point: Vector3<f64>, _normal: Vector3<f64>, _mat: Material, borders: Option<[f64; 6]>, maxt: f64) -> Plane {
        match borders {
            Some(b) => Plane { 
                normal: _normal, 
                holder: vec3_dot(_point, _normal), 
                center: _point, 
                mat: _mat,
                min_x: b[0],
                max_x: b[1],
                min_y: b[2],
                max_y: b[3],
                min_z: b[4],
                max_z: b[5]
            },
            None => Plane { 
                normal: _normal, 
                holder: vec3_dot(_point, _normal), 
                center: _point, 
                mat: _mat,
                min_x: -maxt,
                max_x: maxt,
                min_y: -maxt,
                max_y: maxt,
                min_z: -maxt,
                max_z: maxt
            }
        }
        
    }

    pub fn hit(self, ray: Ray, mint: f64, maxt: f64, rec: &mut HitRecord) -> bool {
        let denominator = vec3_dot(self.normal, ray.direction());

        if denominator != 0.0 {
            let t = (self.holder - vec3_dot(self.normal, ray.origin())) / denominator;
            let p = ray.point_by_t(t);
            if t < maxt && t > mint && p[0] > self.min_x && p[0] < self.max_x && p[1] > self.min_y && p[1] < self.max_y && p[2] > self.min_z && p[2] < self.max_z {
                rec.t = t;
                rec.p = p;
                rec.normal = vec3_normalized(self.normal);
                rec.mat = self.mat;
                rec.center = self.center;
                return true;
            }
        }
        false
    }
}

#[derive(Debug)]
struct HitRecord {
    t: f64,
    p: Vector3<f64>,
    normal: Vector3<f64>,
    mat: Material,
    center: Vector3<f64>
}

impl HitRecord {
    pub fn new() -> HitRecord {
        HitRecord {
            t: 0.0,
            p: [0.0,0.0,0.0],
            normal: [0.0,0.0,0.0],
            mat: Material::standard(),
            center: [0.0,0.0,0.0]
        }
    }
}

#[derive(Clone,Copy,Debug)]
struct Material {
    albedo: Vector3<f64>,
    fuzziness: f64,
    t: MaterialType,
    texture: Option<[[u8;3];2500]>,
    scale: f64
}

impl Material {
    pub fn standard() -> Material{
        Material {
            albedo: [0.0,0.0,0.0],
            fuzziness: 0.0,
            t: MaterialType::Lambertian,
            texture: None,
            scale: 0.0
        }
    }

    pub fn new_l(_albedo: Vector3<f64>) -> Material {
        Material {
            albedo: _albedo,
            fuzziness: 0.0,
            t: MaterialType::Lambertian,
            texture: None,
            scale: 0.0
        }
    }

    pub fn new_m(_albedo: Vector3<f64>, _f: f64) -> Material {
        Material {
            albedo: _albedo,
            fuzziness: _f,
            t: MaterialType::Metal,
            texture: None,
            scale: 0.0
        }
    }

    pub fn new_t(_texture: &str, _scale: f64) -> Material {
        let decoder = png::Decoder::new(File::open(_texture).unwrap());
        let mut reader = decoder.read_info().unwrap();

        let mut buf = vec![0; reader.output_buffer_size()];
        
        let info = reader.next_frame(&mut buf).unwrap();
        let bytes = &buf[..info.buffer_size()];

        let mut arr = [[0; 3]; 2500];
        let step: usize = match bytes.len() {7500 => 3, _ => 4};
        for i in (0..bytes.len()).step_by(step) {
            arr[i/step] = [bytes[i],bytes[i+1],bytes[i+2]];
        }

        Material { 
            albedo: [0.0,0.0,0.0], 
            fuzziness: 0.0, 
            t: MaterialType::Texture,
            texture: Some(arr),
            scale: _scale
        }
    }

    pub fn scatter(self, ray: Ray, rec: HitRecord, attenuation: &mut Vector3<f64>, scattered: &mut Ray) -> bool {
        match self.t {
            MaterialType::Lambertian => {
                *scattered = Ray::new(rec.p, vec3_add(rec.normal, get_random_in_sphere()));
                *attenuation = self.albedo;
                true
            },
            MaterialType::Metal => {
                let mut reflected = vec3_normalized(vec3_sub(ray.direction(),vec3_scale(rec.normal, vec3_dot(rec.normal,ray.direction())*2.0)));
                reflected = vec3_add(reflected,vec3_scale(get_random_in_sphere(), self.fuzziness));
                *scattered = Ray::new(rec.p, reflected);
                *attenuation = self.albedo;
                vec3_dot(scattered.direction(), rec.normal) > 0.0
            },
            MaterialType::Texture => {
                *scattered = Ray::new(rec.p, vec3_add(rec.normal, get_random_in_sphere()));

                let offset = vec3_sub(rec.p, rec.center);
                let horizontal_vector: Vector3<f64>;
                let vertical_vector: Vector3<f64>;
                let stdvecs = [[1.0,0.0,0.0],[-1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,1.0],[0.0,0.0,-1.0]];
                if stdvecs.iter().any(|v| *v == rec.normal) {
                    let (hv, vv) = get_vecs(rec.normal);
                    horizontal_vector = hv;
                    vertical_vector= vv;
                }
                else {
                    horizontal_vector = [1.0,0.0,-(rec.normal[0]/rec.normal[2])];
                    vertical_vector = vec3_cross(rec.normal, horizontal_vector);
                }
                let horizontal_offset = vec3_dot(horizontal_vector, offset) / vec3_len(horizontal_vector);
                let vertical_offset = vec3_dot(vertical_vector, offset) / vec3_len(vertical_vector);
                let x = horizontal_offset / self.scale;
                let y = vertical_offset / self.scale;
                let xs: usize;
                let ys: usize;
                if x < 0.0 {
                    xs = 49-(x.abs() as usize % 50)
                }
                else {
                    xs = x as usize % 50
                }
                if y < 0.0 {
                    ys = 49-(y.abs() as usize % 50)
                }
                else {
                    ys = y as usize % 50
                }
                let pixel = ys * 50 + xs;
                *attenuation = [self.texture.unwrap()[pixel][0] as f64, self.texture.unwrap()[pixel][1] as f64, self.texture.unwrap()[pixel][2] as f64];
                *attenuation = vec3_scale(*attenuation, 1.0/255.0);
                true
            }
        }
    }
}

#[derive(Clone,Copy,Debug)]
enum MaterialType {
    Lambertian,
    Metal,
    Texture
}

fn main() {
    let path = Path::new(r"troll_render.png");
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let width:usize = 600;
    let height: usize = 300;

    /*
    // test
    let spheres = vec![
        Sphere::new([0.0,0.0,-2.0], 1.0, Material::new_m([0.3,0.8,0.3], 0.2)),
        Sphere::new([2.0,-0.5,-1.5], 0.5, Material::new_l([0.1,0.3,0.6])),
        Sphere::new([4.0,0.5,-4.0], 1.5, Material::new_m([0.8,0.8,0.8], 0.05)),
        Sphere::new([4.0,3.0,-4.0], 1.0, Material::new_m([0.3,0.6,0.2], 1.0))
    ];

    let planes = vec![
        Plane::new([0.0,-1.0,-2.0], [0.0,1.0,0.0], Material::new_l([0.8,0.1,0.8])),
        Plane::new([-15.0,0.0,-15.0], [15.0,0.0,10.0], Material::new_t("brick.png", 0.3))
    ];
    */

    // troll
    let spheres = vec![
        Sphere::new([0.0,-22.0,-25.0], 3.0, Material::new_m([0.3,0.8,0.5], 0.6)),
        Sphere::new([6.0,-12.0,-45.0], 3.0, Material::new_m([0.1,0.2,0.5], 0.05)),
        Sphere::new([-22.0,-22.0,-25.0], 3.0, Material::new_m([0.3,0.8,0.5], 0.6)),
        Sphere::new([6.0,-20.0,-45.0], 5.0, Material::new_m([0.9,0.2,0.5], 0.8))
    ];

    let planes = vec![
        Plane::new([-25.0,0.0,0.0], [1.0,0.0,0.0], Material::new_t("brick.png", 0.5), Some([-1000.0, 1000.0, -25.0, 25.0, -50.0, 0.0]), 1000.0),
        Plane::new([25.0,25.0,0.0], [-1.0,0.0,0.0], Material::new_t("troll.png", 1.0), Some([-1000.0, 1000.0, -25.0, 25.0, -50.0, 0.0]), 1000.0),
        Plane::new([-25.0,25.0,-50.0], [0.0,0.0,1.0], Material::new_m([1.0,1.0,1.0], 0.0), Some([-25.0,25.0,-25.0,25.0,-1000.0,1000.0]), 1000.0),
        Plane::new([-25.0,25.0,0.0], [0.0,0.0,-1.0], Material::new_m([1.0,1.0,1.0], 0.0), Some([-25.0,25.0,-25.0,25.0,-1000.0,1000.0]), 1000.0),
        Plane::new([0.0,-25.0,0.0], [0.0,1.0,0.0], Material::new_l([1.0,0.2,0.8]), None, 1000.0)
    ];

    let maxt = 1000.0;
    let mint = 0.001;

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
    let mut pixeles = create_rays(Camera::new([0.0,0.0,0.0]), width, height);
    trace_all(spheres, planes,&mut pixeles,&mut data, mint, maxt);

    
    writer.write_image_data(&data).unwrap();
}

fn create_rays(_camera: Camera, _w: usize, _h: usize) -> Vec<[Ray;25]>{
    let mut rays: Vec<[Ray;25]> = Vec::new();
    let mut rng = rand::thread_rng();

    let top_left = vec3_sub(_camera.tl,_camera.c);
    let width_offset = vec3_sub(_camera.tr,_camera.tl);
    let height_offset = vec3_sub(_camera.bl,_camera.tl);

    for i in 0.._h {
        for j in 0.._w {
            let mut n: Vec<f64> = Vec::new();
            for _ in 0..50 {
                n.push(rng.gen());
            }
            let mut vec: Vec<Ray> = Vec::new();
            for k in 0..25 {
                vec.push(Ray::new(_camera.c, vec3_add(top_left, vec3_add(vec3_scale(width_offset, (j as f64 + n[k*2])/_w as f64), vec3_scale(height_offset, (i as f64 + n[k*2+1])/_h as f64)))));
            }
            let mut arr: [Ray; 25] = [Ray::new([0.0,0.0,0.0], [0.0,0.0,0.0]); 25];
            for k in 0..25 {
                arr[k] = vec[k];
            }
            rays.push(arr);
        }
    }

    rays
}


fn trace_all(spheres: Vec<Sphere>, planes: Vec<Plane>, pixeles: &mut Vec<[Ray;25]>, data: &mut [u8; 600 * 300 * 4], mint: f64, maxt: f64) {
    
    let pix_c: Vec<[f64; 3]> = (0..pixeles.len()).into_par_iter().map(|i| {
        let mut pixel_colour: [f64; 3] = [0.0,0.0,0.0];
        for j in pixeles[i] {
            pixel_colour = vec3_add(pixel_colour, trace(j, &spheres, &planes, mint, maxt, 0));
        }
        let mut pc  = vec3_scale(pixel_colour, 1.0/25.0);
        pc = [pc[0].sqrt(), pc[1].sqrt(),pc[2].sqrt()];
        vec3_scale(pc, 255.0)
    }).collect();
    
    for i in 0..pix_c.len() {
        
        assign_pixel(i, pix_c[i][0] as u8, pix_c[i][1] as u8, pix_c[i][2] as u8, 255, data);
    }
    
    /*
    for i in 0..pixeles.len() {
        let mut pixel_colour: [f64; 3] = [0.0,0.0,0.0];
        for j in pixeles[i] {
            pixel_colour = vec3_add(pixel_colour, trace(j, &spheres, &planes, mint, maxt, 0));
        }
        let mut pc  = vec3_scale(pixel_colour, 1.0/25.0);
        pc = [pc[0].sqrt(), pc[1].sqrt(),pc[2].sqrt()];
        pc = vec3_scale(pc, 255.0);
        assign_pixel(i, pc[0] as u8, pc[1] as u8, pc[2] as u8, 255, data);
    }
    */
}

fn trace(ray: Ray, spheres: &Vec<Sphere>, planes: &Vec<Plane>, mint: f64, maxt: f64, depth: u8) -> Vector3<f64> {
    let mut rec: HitRecord = HitRecord::new();
    if hits(ray, spheres, planes, mint, maxt, &mut rec) {
        let mut scattered: Ray = Ray::new([0.0,0.0,0.0], [0.0,0.0,0.0]);
        let mut attenuation: Vector3<f64> = [0.0,0.0,0.0];
        if depth < 50 && rec.mat.scatter(ray,rec,&mut attenuation,&mut scattered) { 
            return vec3_mul(trace(scattered,spheres, planes,mint,maxt, depth + 1), attenuation);
        }
        else {
            return [0.0,0.0,0.0];
        }
    }
    else {
        let t = 0.5 * (vec3_normalized(ray.direction())[1] + 1.0);
        return vec3_add(vec3_scale([1.0,1.0,1.0], 1.0-t),vec3_scale([0.2,0.4,1.0], t));
    }
}

fn hits(ray: Ray, spheres: &Vec<Sphere>, planes: &Vec<Plane>, mint: f64, maxt: f64, rec: &mut HitRecord) -> bool {
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
            rec.mat = temp_rec.mat;
        }
    }
    for k in 0..planes.len(){
        if planes[k].hit(ray,mint,min_dist,&mut temp_rec) {
            hit_anything = true;
            min_dist = temp_rec.t;
            rec.t = temp_rec.t;
            rec.p = temp_rec.p;
            rec.normal = temp_rec.normal;
            rec.mat = temp_rec.mat;
            rec.center = temp_rec.center;
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

fn get_vecs(normal: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    let empty = normal.iter().position(|c| *c != 0.0).unwrap();
    let mut vecs: [Vector3<f64>; 2] = [[0.0,0.0,0.0]; 2];
    let mut j: usize = 0;
    for i in [0,2,1] {
        if i != empty {
            vecs[j][i] = 1.0;
            j += 1;
        }
    }
    (vecs[0],vecs[1])
}

fn assign_pixel(pixel: usize, r: u8, g: u8, b: u8, a: u8, data: &mut [u8; 600 * 300 * 4]) {
    let i = pixel * 4;
    data[i] = r;
    data[i+1] = g;
    data[i+2] = b;
    data[i+3] = a;
}
