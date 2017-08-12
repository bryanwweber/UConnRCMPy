"""
Test module for the simulations module
"""
import numpy as np
import pytest
from ..simulations import Simulation, VolumeProfile


class TestVolumeProfile(object):
    """Tests for the VolumeProfile"""
    @pytest.fixture
    def vp(self):
        time = np.array([0, 1, 2, 3, 4])
        volume = np.array([10, 20, 40, 70, 110])
        return VolumeProfile(time, volume)

    def test_volume_normalized(self, vp):
        """Test that the volume is normalized"""
        assert np.all(vp.volume == np.array([1.0, 2.0, 4.0, 7.0, 11.0]))

    def test_velocity_correct(self, vp):
        """Test the velocity calculation is correct"""
        assert np.all(vp.velocity == np.array([1.0, 2.0, 3.0, 4.0, 0.0]))

    def test_vel_at_time_gt_end(self, vp):
        """Test that velocities after the end time are zero"""
        assert vp(5) == 0.0

    def test_vel_at_int_time(self, vp):
        assert vp(1.0) == 2.0
        assert vp(1.5) == 2.0
        assert vp(2.0) == 3.0
        assert vp(2.5) == 3.0
