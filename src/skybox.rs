use bevy::{
    core_pipeline::Skybox,
    prelude::*,
    render::render_resource::{Extent3d, TextureViewDescriptor, TextureViewDimension},
};
use image::{EncodableLayout, ImageBuffer, Rgba32FImage};

const IMAGE_RES: u64 = 3700;

/// Plugin for the Star Skybox
pub struct StarSkyboxPlugin;

impl Plugin for StarSkyboxPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(SkyboxImage(None))
            .add_systems(Startup, generate_skybox)
            .add_systems(Update, attach_skybox_image_to_cameras);
    }
}

/// Marker component to indicate which camera should have the skybox background
#[derive(Debug, Component)]
pub struct StarSkybox;

/// Resource to hold the generated image of the skybox
#[derive(Debug, Default, Resource)]
pub struct SkyboxImage(pub Option<Handle<Image>>);

fn generate_skybox(mut skybox_image_res: ResMut<SkyboxImage>, mut images: ResMut<Assets<Image>>) {
    let image = generate_skybox_image();
    let mut image_asset: Image = Image::new(
        Extent3d {
            width: image.width(),
            height: image.height(),
            depth_or_array_layers: 1,
        },
        bevy::render::render_resource::TextureDimension::D2,
        image.as_bytes().to_vec(),
        bevy::render::render_resource::TextureFormat::Rgba32Float,
    );
    image_asset.reinterpret_stacked_2d_as_array(
        image_asset.texture_descriptor.size.height / image_asset.texture_descriptor.size.width,
    );
    image_asset.texture_view_descriptor = Some(TextureViewDescriptor {
        dimension: Some(TextureViewDimension::Cube),
        ..default()
    });
    let handle = images.add(image_asset);
    skybox_image_res.0 = Some(handle);
}

fn attach_skybox_image_to_cameras(
    mut commands: Commands,
    mut cameras_query: Query<Entity, (With<StarSkybox>, Without<Skybox>)>,
    skybox_image_res: Res<SkyboxImage>,
) {
    // If the skybox image isn't generated yet, don't do anything
    if skybox_image_res.0.is_none() {
        return;
    }

    // Attach the skybox image to all cameras with the StarSkybox component and no skybox
    for entity in cameras_query.iter_mut() {
        commands
            .entity(entity)
            .insert(Skybox(skybox_image_res.0.clone().unwrap()));
    }
}

fn generate_skybox_image() -> Rgba32FImage {
    let stars = process_hipparcos_catalog(&hipparcos_catalog::STARS);

    // pos x
    let pos_x_stars = stars
        .iter()
        .filter(|star| matches!(star.face, Face::PositiveX))
        .collect::<Vec<&Star>>();
    let pos_x_image = stars_to_image(pos_x_stars, IMAGE_RES);

    // neg x
    let neg_x_stars = stars
        .iter()
        .filter(|star| matches!(star.face, Face::NegativeX))
        .collect::<Vec<&Star>>();
    let neg_x_image = stars_to_image(neg_x_stars, IMAGE_RES);

    // pos y
    let pos_y_stars = stars
        .iter()
        .filter(|star| matches!(star.face, Face::PositiveY))
        .collect::<Vec<&Star>>();
    let pos_y_image = stars_to_image(pos_y_stars, IMAGE_RES);

    // neg y
    let neg_y_stars = stars
        .iter()
        .filter(|star| matches!(star.face, Face::NegativeY))
        .collect::<Vec<&Star>>();
    let neg_y_image = stars_to_image(neg_y_stars, IMAGE_RES);

    // pos z
    let pos_z_stars = stars
        .iter()
        .filter(|star| matches!(star.face, Face::PositiveZ))
        .collect::<Vec<&Star>>();
    let pos_z_image = stars_to_image(pos_z_stars, IMAGE_RES);

    // neg z
    let neg_z_stars = stars
        .iter()
        .filter(|star| matches!(star.face, Face::NegativeZ))
        .collect::<Vec<&Star>>();
    let neg_z_image = stars_to_image(neg_z_stars, IMAGE_RES);

    // Concatenate the images into a single image
    // pos_x_image
    // neg_x_image
    // pos_y_image
    // neg_y_image
    // neg_z_image
    // pos_z_image
    let mut buf = ImageBuffer::new(IMAGE_RES as u32, (IMAGE_RES * 6) as u32);

    for (x, y, p) in buf.enumerate_pixels_mut() {
        let x = x as u64;
        let y = y as u64;

        let mut pixel = image::Rgba::<f32>([0.0, 0.0, 0.0, 1.0]);

        if y < IMAGE_RES {
            // pos x
            pixel = *pos_x_image.get_pixel(x as u32, y as u32);
        } else if y < IMAGE_RES * 2 {
            // neg x
            pixel = *neg_x_image.get_pixel(x as u32, (y - IMAGE_RES) as u32);
        } else if y < IMAGE_RES * 3 {
            // pos y
            pixel = *pos_y_image.get_pixel(x as u32, (y - IMAGE_RES * 2) as u32);
        } else if y < IMAGE_RES * 4 {
            // neg y
            pixel = *neg_y_image.get_pixel(x as u32, (y - IMAGE_RES * 3) as u32);
        } else if y < IMAGE_RES * 5 {
            // neg z
            pixel = *neg_z_image.get_pixel(x as u32, (y - IMAGE_RES * 4) as u32);
        } else if y < IMAGE_RES * 6 {
            // pos z
            pixel = *pos_z_image.get_pixel(x as u32, (y - IMAGE_RES * 5) as u32);
        }

        *p = pixel
    }

    buf
}

#[derive(Debug)]
struct Star {
    face: Face,
    x: u64,
    y: u64,
    luminance: f32,
}

fn process_hipparcos_catalog(star_catalog: &[hipparcos_catalog::Star]) -> Vec<Star> {
    let mut parsed_stars = vec![];

    for star in star_catalog {
        let mag = star.magnitude;
        let ra = star.ra_deg;
        let dec = star.dec_deg;

        // Calculate the luminance from the magnitude (https://en.wikipedia.org/wiki/Apparent_magnitude#Calculating_a_star's_apparent_magnitude)
        let luminance: f32 = 10f32.powf(-0.1 * mag as f32);

        // Convert the RA/DEC to a cubemap face and x-y coordinates
        let (face, x, y) = ra_dec_to_face(ra, dec, IMAGE_RES);

        let star = Star {
            face,
            x: x as u64,
            y: y as u64,
            luminance,
        };

        parsed_stars.push(star);
    }

    parsed_stars
}

#[derive(Debug)]
enum Face {
    PositiveX,
    NegativeX,
    PositiveY,
    NegativeY,
    PositiveZ,
    NegativeZ,
}

fn ra_dec_to_face(ra: f64, dec: f64, res: u64) -> (Face, f64, f64) {
    // Convert the degs RA/DEC to a cubemap face and x-y coordinates
    // Convert RA/DEC to Cartesian Coordinates
    let ra_rad = ra * std::f64::consts::PI / 180.0;
    let dec_rad = dec * std::f64::consts::PI / 180.0;

    let x = dec_rad.cos() * ra_rad.cos();
    let y = dec_rad.sin();
    let z = -dec_rad.cos() * ra_rad.sin();

    // Determine the Cube Face: For each Cartesian coordinate (X, Y, Z), find which value has the maximum absolute value. The sign of this value will indicate which face of the cube the star should be projected onto:
    let mut face;
    // If abs(X) is the greatest, then the star projects onto the left face if X is negative and the right face if X is positive.
    let mut max = x.abs();
    if x < 0.0 {
        face = Face::NegativeX;
    } else {
        face = Face::PositiveX;
    }

    // If abs(Y) is the greatest, then the star projects onto the down face if Y is negative and the up face if Y is positive.
    if y.abs() > max {
        max = y.abs();
        if y < 0.0 {
            face = Face::NegativeY;
        } else {
            face = Face::PositiveY;
        }
    }
    // If abs(Z) is the greatest, then the star projects onto the back face if Z is negative and the front face if Z is positive.
    if z.abs() > max {
        if z < 0.0 {
            face = Face::NegativeZ;
        } else {
            face = Face::PositiveZ;
        }
    }

    // Map the Coordinates to a 2D Point on the Cube Face: Normalize the remaining two coordinates (those other than the one with the maximum absolute value) to be within [-1, 1], then map these normalized coordinates to the pixel resolution of your cubemap image.

    // For instance, if the star projects onto the right face (positive X), calculate the normalized 2D coordinates as follows:
    // U = Y / abs(X) / 2 + 0.5
    // V = Z / abs(X) / 2 + 0.5
    let u = match face {
        Face::PositiveX | Face::NegativeX => y / x.abs() / 2.0 + 0.5,
        Face::PositiveY | Face::NegativeY => x / y.abs() / 2.0 + 0.5,
        Face::PositiveZ | Face::NegativeZ => x / z.abs() / 2.0 + 0.5,
    };
    let v = match face {
        Face::PositiveX | Face::NegativeX => z / x.abs() / 2.0 + 0.5,
        Face::PositiveY | Face::NegativeY => z / y.abs() / 2.0 + 0.5,
        Face::PositiveZ | Face::NegativeZ => y / z.abs() / 2.0 + 0.5,
    };

    // Map the 2D Point to a Pixel in the Cubemap Image: Multiply the normalized coordinates by the resolution of the cubemap image to get the pixel coordinates of the star in the cubemap image.
    let x = u * res as f64;
    let y = (1.0 - v) * res as f64;

    (face, x, y)
}

fn stars_to_image(stars: Vec<&Star>, res: u64) -> Rgba32FImage {
    let mut image = Rgba32FImage::new(res as u32, res as u32);

    // Fill with black
    for (_, _, p) in image.enumerate_pixels_mut() {
        *p = image::Rgba::<f32>([0.0, 0.0, 0.0, 1.0]);
    }

    for star in stars {
        let x = star.x as u32;
        let y = star.y as u32;
        let luminance = star.luminance;

        image.put_pixel(
            x,
            y,
            image::Rgba::<f32>([luminance, luminance, luminance, 1.0]),
        );
    }

    image
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_skybox() {
        let skybox = generate_skybox_image();
        assert_eq!(skybox.width(), IMAGE_RES as u32);
        assert_eq!(skybox.height(), IMAGE_RES as u32 * 6);
    }
}
