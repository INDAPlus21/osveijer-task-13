use std::path::Path;
use std::fs::File;
use std::io::{BufWriter};
use png;

fn main() {
    let path = Path::new(r"test.png");
    let file = File::create(path).unwrap();
    let ref mut w = BufWriter::new(file);

    let width:usize = 200;
    let height: usize = 200;

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

    let mut data: [u8; 200 * 200 * 4] = [0; 200 * 200 * 4];
    for i in 0..200 {
        for j in 0..200 {
            assign_pixel(i*200+j, ((i as f32/200.0)*255.0) as u8, ((j as f32/200.0)*255.0) as u8, 0, 255, &mut data)
        }
    }
    writer.write_image_data(&data).unwrap();
}

fn assign_pixel(pixel: usize, r: u8, g: u8, b: u8, a: u8, data: &mut [u8; 200 * 200 * 4]) {
    let i = pixel * 4;
    data[i] = r;
    data[i+1] = g;
    data[i+2] = b;
    data[i+3] = a;
}
