#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include <hip/hip_runtime.h>

#include "definitions.h"
#include "config.h"
#include "particle.h"
#include "beam_elements.h"
#include "beamfields.h"
#include "lattice.h"

__global__ 
//__attribute__((amdgpu_flat_work_group_size(DEMOTRACK_DEFAULT_BLOCK_SIZE, DEMOTRACK_DEFAULT_BLOCK_SIZE)))
__launch_bounds__(DEMOTRACK_DEFAULT_BLOCK_SIZE, 3)	
//__attribute__((amdgpu_waves_per_eu(DEMOTRACK_DEFAULT_BLOCK_SIZE/64)))
void Track_particles_until_turn01_01(
    demotrack::Particle* particle_set,
    demotrack::int64_type const num_particles,
    double const* __restrict__ lattice_buffer,
    demotrack::uint64_type const max_lattice_buffer_index,
    demotrack::int64_type const until_turn )
{
    namespace dt = demotrack;

    dt::int64_type const STRIDE = blockDim.x * gridDim.x;
    dt::int64_type idx = threadIdx.x + blockIdx.x * blockDim.x;

    if( idx == 0 ) printf( "info :: beam-fields enabled in kernel\r\n" );

    for( ; idx < num_particles ; idx += STRIDE )
    {
        dt::Particle p = particle_set[ idx ];
        dt::uint64_type const start_at_element = p.at_element;

        while( ( p.state == 1 ) && ( p.at_turn < until_turn ) )
        {
            dt::uint64_type slot_idx = 0;

            while( ( p.state == 1 ) && ( slot_idx < max_lattice_buffer_index ) )
            {
                /* all elements are stored with their type_id as the first
                 * data member -> retrieve this number and dispatch
                 * the track method accordingly */

                dt::beam_element_type const type_id = ( dt::beam_element_type )(
                    int )lattice_buffer[ slot_idx ];

                switch( type_id )
                {
                    case dt::BEAM_ELEMENT_DRIFT: // cf. beam_elements.h
                    {
                        const dt::Drift *const __restrict__ elem =
                            ( dt::Drift const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        dt::Drift::GLOBAL_APERTURE_CHECK( p );
                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_DRIFT_EXACT: // cf. beam_elements.h
                    {
                        const dt::DriftExact *const __restrict__ elem =
                            ( dt::DriftExact const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        // Use GLOBAL_APERTURE_CHECK from Drift -> it's the same
                        dt::Drift::GLOBAL_APERTURE_CHECK( p );
                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_MULTIPOLE: // cf. beam_elements.h
                    {
                        const dt::Multipole *const __restrict__ elem =
                            ( dt::Multipole const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_XY_SHIFT: // cf. beam_elements.h
                    {
                        const dt::XYShift *const __restrict__ elem =
                            ( dt::XYShift const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_SROTATION: // cf. beam_elements.h
                    {
                        const dt::SRotation *const __restrict__ elem =
                            ( dt::SRotation const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_CAVITY: // cf. beam_elements.h
                    {
                        const dt::Cavity *const __restrict__ elem =
                            ( dt::Cavity const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_LIMIT_RECT: // cf. beam_elements.h
                    {
                        const dt::LimitRect *const __restrict__ elem =
                            ( dt::LimitRect const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_LIMIT_ELLIPSE: // cf. beam_elements.h
                    {
                        const dt::LimitEllipse *const __restrict__ elem =
                            ( dt::LimitEllipse const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_LIMIT_RECT_ELLIPSE: // cf. beam_elements.h
                    {
                        const dt::LimitRectEllipse *const __restrict__ elem =
                            ( dt::LimitRectEllipse const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_DIPEDGE: // cf. beam_elements.h
                    {
                        const dt::DipoleEdge *const __restrict__ elem =
                            ( dt::DipoleEdge const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    case dt::BEAM_ELEMENT_SC_COASTING: // cf. beamfields.h
                    {
                        const dt::SpaceChargeCoasting *const __restrict__ elem =
                            ( dt::SpaceChargeCoasting const* )&lattice_buffer[ slot_idx ];

                        dt::uint64_type const next_slot_idx =
                            elem->track( p, slot_idx );

                        slot_idx = next_slot_idx;
                        break;
                    }

                    default:
                    {
                        /* unknown beam element -> loose particle and quit */
                        p.state = 0;
                        slot_idx = max_lattice_buffer_index;
                    }
                };

            }

            if( p.state == 1 )
            {
                p.at_element = start_at_element;
                ++p.at_turn;
            }
        }

        particle_set[ idx ] = p;
    }
}

