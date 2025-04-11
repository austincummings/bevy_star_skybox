use bevy::{
    core_pipeline::{
        auto_exposure::{AutoExposure, AutoExposurePlugin},
        bloom::Bloom,
    },
    prelude::*,
};
use bevy_star_skybox::{StarSkybox, StarSkyboxPlugin};

fn main() {
    App::new()
        .add_plugins((DefaultPlugins, AutoExposurePlugin, StarSkyboxPlugin))
        .add_systems(Startup, setup)
        .run();
}

fn setup(mut commands: Commands) {
    commands.spawn((
        Transform::from_xyz(0.0, 0.0, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        Camera3d::default(),
        Camera {
            hdr: true,
            ..default()
        },
        Bloom::default(),
        AutoExposure::default(),
        StarSkybox,
    ));
}
