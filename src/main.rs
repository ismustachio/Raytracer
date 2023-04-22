extern crate math_engine;
use std::fs::File;
use std::io::Write;
use math_engine::prelude::*;
use std::{thread, time};
use rand::Rng;
use rand::rngs::ThreadRng;
use rand::distributions::{Distribution, Uniform};

pub struct Camera {
    pub origin: Point3,
    pub lower_left_corner_point: Point3,
    pub vertical: Vector3,
    pub horizontal: Vector3,
}

impl Camera {
    fn new() -> Camera {
        let aspect_ratio = 16.0 / 9.0;
        let viewport_height = 2.0;
        let viewport_width = aspect_ratio * viewport_height;
        let focal_length = 1.0; 
        let origin = Default::default();
        let horizontal = Vector3::new(viewport_width, 0.0, 0.0);
        let vertical = Vector3::new(0.0, viewport_height, 0.0);
        let lower_left_corner_point = origin - (horizontal / 2.0) - (vertical / 2.0) - Vector3::new(0.0, 0.0, focal_length);
        Camera {
            origin,
            lower_left_corner_point,
            vertical,
            horizontal,
        }
    }

    fn get_ray(&self, u: f32, v: f32) -> Ray {
        Ray::new(self.origin, self.lower_left_corner_point + (self.horizontal * u) + (self.vertical*v) - self.origin)
    }
}

fn random_vector3(rng: &mut ThreadRng) -> Vector3 {
    let roll = Uniform::from(0.0..1.0); 
    Vector3::new(roll.sample(rng), roll.sample(rng), roll.sample(rng))
}

fn random_vector3_range(rng: &mut ThreadRng, min: f32, max: f32) -> Vector3 {
    let roll = Uniform::from(min..max); 
    Vector3::new(roll.sample(rng), roll.sample(rng), roll.sample(rng))
}

fn random_normalized_vector3(rng: &mut ThreadRng) -> Vector3 {
    loop {
        let v = random_vector3_range(rng, -1.0, 1.0); 
        if v.dot(&v) >= 1.0 {
            continue;
        }
        return v.normalize();
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Sphere {
    pub id: u64,
    pub center: Point3,
    pub radius: f32, 
}

impl Sphere {
    fn new(center: Point3, radius: f32) -> Sphere {
        let id = 1;
        Sphere {
            id,
            center,
            radius,
        }
    }
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32, hit_record: &mut HitRecord) -> bool {
        let sphere_to_ray = ray.origin - self.center;
        let a = ray.direction.dot(&ray.direction);
        let half_b = sphere_to_ray.dot(&ray.direction);
        let c = sphere_to_ray.dot(&sphere_to_ray) - ( self.radius * self.radius );
        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 { return false; }
        let sqrtd = discriminant.sqrt();
        let mut root = (-half_b - sqrtd) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrtd) / a;
            if root < t_min || t_max < root {
                return false;
            }
        }
        hit_record.t = root;
        hit_record.position = ray.position(root);
        let normal = (hit_record.position - self.center) / self.radius;
        hit_record.set_face_normal(ray, normal);
        return true;
    }

    fn id(&self) -> u64 {
        self.id
    }
}

struct Hittables<T: Hittable> {
    pub hittables: Vec<T>,
    id: u64,
}

impl<T: Hittable> Hittables<T> {
    fn new(hittables: Vec<T>) -> Hittables<T> {
        let id = 1;
        Hittables {
            hittables,
            id,
        }
    }
}

impl<T: Hittable> Hittable for Hittables<T> {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32, hit_record: &mut HitRecord) -> bool {
        let mut result = false;
        let mut closest_so_far = t_max;
        let mut temp_record = Default::default();
        for hittable in &self.hittables {
            if hittable.hit(ray, t_min, closest_so_far, &mut temp_record) {
               result = true;
               closest_so_far = temp_record.t;
               *hit_record = temp_record.clone();
            }
        }
        result
    }

    fn id(&self) -> u64 {
        self.id
    }
}

#[derive(Debug, Default, Copy, Clone)]
struct HitRecord {
    pub position: Point3,
    pub normal: Vector3,
    pub t: f32,
    pub front_face: bool,
}

impl HitRecord {
    fn new(position: Point3, normal: Vector3, t: f32, front_face: bool) -> HitRecord {
        HitRecord{
            position,
            normal,
            t,
            front_face,
        }
    }

    fn set_face_normal(&mut self, ray: &Ray, normal: Vector3) {
        self.front_face = ray.direction.dot(&normal) < 0.0;
        if self.front_face {
            self.normal = normal;
        } else {
            self.normal = Vector3::new(-normal.x, -normal.y, -normal.z);
        }
    }
}

trait Hittable {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32, hit_record: &mut HitRecord) -> bool;

    fn id(&self) -> u64;
}


#[derive(Debug, Default, Copy, Clone)]
pub struct Ray {
    pub origin: Point3,
    pub direction: Vector3,
}

impl Ray {
    fn new(origin: Point3, direction: Vector3) -> Ray {
        Ray{
            origin,
            direction,
        }
    }

    fn position(&self, t: f32) -> Point3 {
        self.origin + (self.direction * t)
    }

    fn color<T: Hittable>(&self, hittable: &T, depth: i32) -> RGB {
        if depth <= 0 {
            return BLACK;
        }

        let mut record = Default::default();
        if hittable.hit(self, 0.001, std::f32::INFINITY, &mut record) {
            let mut rng = rand::thread_rng();
            let target = record.position + record.normal + random_normalized_vector3(&mut rng);
            return Ray::new(record.position, target - record.position).color(hittable, depth - 1) * 0.5;
            //return (RGB::new(record.normal.x, record.normal.y, record.normal.z) + 
            //WHITE) * 0.5;
        }
        let unit_direction = self.direction.normalize();
        let t = 0.5 * (unit_direction.y + 1.0);
        let a = WHITE * (1.0 - t);
        let b = RGB::new(0.5, 0.7, 1.0) * t;
        return a + b;
    }
}

fn clamp(x: f32, min: f32, max: f32) -> f32 {
    if x < min {
        return min;
    }
    if x > max {
        return max;
    }
    return x;
}


fn main() {
    let aspect_ratio = 16.0 / 9.0;
    let max_depth = 50;
    
    // image
    let image_width = 400;
    let image_height = (image_width as f32 / aspect_ratio) as usize;

    // camera
    //let viewport_height = 2.0;
    //let viewport_width = aspect_ratio * viewport_height;
    //let focal_length = 1.0; 
    //let origin = Default::default();
    //let horizontal = Vector3::new(viewport_width, 0.0, 0.0);
    //let vertical = Vector3::new(0.0, viewport_height, 0.0);
    //let mut horizontal_scale = horizontal / 2.0;
    //let mut vertical_scale = vertical / 2.0;
    //let lower_left_corner_point = origin - horizontal_scale - vertical_scale - Vector3::new(0.0, 0.0, focal_length);
    let camera = Camera::new();
    let samples_per_pixel = 100;


    let s1 = Sphere::new(Point3::new(0.0, 0.0, -1.0), 0.5);
    let s2 = Sphere::new(Point3::new(0.0, -100.5, -1.0), 100.0);
    let world = Hittables::new([s1, s2].to_vec());

    let mut rng = rand::thread_rng();
    let roll = Uniform::from(0.0..1.0);

    // render
    let rows = (image_height - 1) as f32;
    let cols = (image_width - 1) as f32;
    let filename = "./t.ppm";
    let ten_millis = time::Duration::from_millis(10);
    let mut buf = format!("P3\n{} {}\n255\n", image_width, image_height);
    for row in (0..=rows as i32).rev() {
        print!("\rScanline remaining: {} ", row);
        std::io::stdout().flush();
        //thread::sleep(ten_millis);
        for col in 0..image_width {
            let mut pixel_color: RGB = Default::default();
            for sample in 0..samples_per_pixel {
                let u = (col as f32 + roll.sample(&mut rng)) / cols;
                let v = (row as f32 + roll.sample(&mut rng)) / rows;
                let ray = camera.get_ray(u, v);
                pixel_color += ray.color(&world, max_depth);
            }
            //let u = row as f32 / rows;
            //let v = col as f32 / cols;
            //let direction = lower_left_corner_point + (horizontal*v) + (vertical*u) - origin;
            //let ray = Ray::new(origin, direction);
            //let mut pixel_color = ray.color(&world);
            // divide the color by the number of samples
            // for anti-aliasing

            // gamma correct for 2.0
            let scale = 1.0 / samples_per_pixel as f32;
            pixel_color.r = (scale * pixel_color.r).sqrt();
            pixel_color.g = (scale * pixel_color.g).sqrt();
            pixel_color.b = (scale * pixel_color.b).sqrt();
            //println!("{:?}", pixel_color);

            buf.push_str(format!("{} {} {}\n", (256.0 * clamp(pixel_color.r, 0.0, 0.999)) as i32,
                                               (256.0 * clamp(pixel_color.g, 0.0, 0.999)) as i32,
                                               (256.0 * clamp(pixel_color.b, 0.0, 0.999)) as i32).as_str());
        }
    }
    buf.push_str("\n");
    println!("\nDone");
    let mut file = File::create(filename).expect("failed to create file");
    file.write_all(buf.as_bytes()).expect("failed to write to file");
}

