use bevy::{core_pipeline::bloom::BloomSettings, prelude::*, render::view::ColorGrading};
use bevy_panorbit_camera::{PanOrbitCamera, PanOrbitCameraPlugin};
use bevy_star_skybox::{StarSkybox, StarSkyboxPlugin};

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, StarSkyboxPlugin))
        .add_plugins(PanOrbitCameraPlugin)
        .add_systems(Startup, setup)
        .run();
}

fn setup(mut commands: Commands) {
    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_xyz(0.0, 0.0, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
            color_grading: ColorGrading {
                exposure: 8.0,
                ..default()
            },
            camera: Camera {
                hdr: true,
                ..default()
            },
            ..default()
        },
        PanOrbitCamera::default(),
        BloomSettings::default(),
        StarSkybox,
    ));
}
