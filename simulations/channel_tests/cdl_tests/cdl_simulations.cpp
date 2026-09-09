#include "cdl_simulations.hpp"

namespace channels::cdl {

    void runSimulation1(void) {

        nrCDLChannel::Config cfg;

        /******************* Predefined Delay Profile ************************/ 
        cfg.DelayProfile.DelayProfile = pdp::PDP::CDL_D;
        cfg.DelayProfile.DelaySpread = 3.0e-8f;

        //   nullopt  -> K-factor scaling disabled
        //   value    -> K-factor scaling enabled
        cfg.DelayProfile.KFactor = std::nullopt;

        //   nullopt         -> angle scaling disabled
        //   MeanAnglesList  -> angle scaling enabled
        cfg.DelayProfile.MeanAngles = std::nullopt;

        /******************* Transmit Antenna Array ************************/ 
        auto& tx = cfg.TransmitAntennaArraySetup;

        // MATLAB: antTx_struct.Size = [2 2 2 1 1];
        tx.ArraySize.M  = 2;
        tx.ArraySize.N  = 2;
        tx.ArraySize.P  = 2;
        tx.ArraySize.Mg = 1;
        tx.ArraySize.Ng = 1;

        // MATLAB: antTx_struct.ElementSpacing = [0.5 0.5 1.0 1.0];
        tx.ElementSpacing.d_v  = 0.5f;
        tx.ElementSpacing.d_h  = 0.5f;
        tx.ElementSpacing.dg_v = 1.0f;
        tx.ElementSpacing.dg_h = 1.0f;

        // MATLAB: antTx_struct.PolarizationAngles = [45 -45];
        tx.PolarizationAngles.theta =  45.0f;
        tx.PolarizationAngles.rho   = -45.0f;

        // MATLAB: antTx_struct.Orientation = [0; 0; 0];
        tx.Orientation.alpha = 0.0f;
        tx.Orientation.beta  = 0.0f;
        tx.Orientation.gamma = 0.0f;

        // MATLAB: antTx_struct.Element = '38.901';
        tx.FieldEffect.Element = antenna::ElementPatternOption::TR_38_901;

        // MATLAB: antTx_struct.PolarizationModel = 'Model-2';
        tx.FieldEffect.PolarizationModel = antenna::PolarizationModelOption::MODEL2;

        /******************* Receive Antenna Array ************************/
        auto& rx = cfg.ReceiveAntennaArraySetup;

        // MATLAB: antRx_struct.Size = [1 1 2 1 1];
        rx.ArraySize.M  = 1;
        rx.ArraySize.N  = 1;
        rx.ArraySize.P  = 2;
        rx.ArraySize.Mg = 1;
        rx.ArraySize.Ng = 1;

        // MATLAB: antRx_struct.ElementSpacing = [0.5 0.5 0.5 0.5];
        rx.ElementSpacing.d_v  = 0.5f;
        rx.ElementSpacing.d_h  = 0.5f;
        rx.ElementSpacing.dg_v = 0.5f;
        rx.ElementSpacing.dg_h = 0.5f;

        // MATLAB: antRx_struct.PolarizationAngles = [0 90];
        rx.PolarizationAngles.theta = 0.0f;
        rx.PolarizationAngles.rho   = 90.0f;

        // MATLAB: antRx_struct.Orientation = [0; 0; 0];
        rx.Orientation.alpha = 0.0f;
        rx.Orientation.beta  = 0.0f;
        rx.Orientation.gamma = 0.0f;

        // MATLAB: antRx_struct.Element = 'isotropic';
        rx.FieldEffect.Element = antenna::ElementPatternOption::ISOTROPIC;

        // MATLAB: antRx_struct.PolarizationModel = 'Model-2';
        rx.FieldEffect.PolarizationModel = antenna::PolarizationModelOption::MODEL2;

        /******************* Construct Channel Object ************************/ 
        //nrCDLChannel cdl;
        //cdl.configure(cfg); 
        // OR
        nrCDLChannel cdl{cfg};


    }
    
}