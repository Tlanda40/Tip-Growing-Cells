import bpy
import bmesh
from mathutils import Vector

# --- 1. Clean up scene ---
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# --- 2. Add UV Sphere (for hemisphere) ---
bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=16, radius=1, location=(0, 0, 0))
hemi_obj = bpy.context.active_object
hemi_mesh = hemi_obj.data

# Enter Edit Mode to delete bottom half
bpy.ops.object.mode_set(mode='EDIT')
bm = bmesh.from_edit_mesh(hemi_mesh)

# Delete vertices below z=0
for v in bm.verts:
    v.select = v.co.z < 0
bmesh.ops.delete(bm, geom=[v for v in bm.verts if v.select], context='VERTS')
bmesh.update_edit_mesh(hemi_mesh)
bpy.ops.object.mode_set(mode='OBJECT')

# --- 3. Add Circle (flat disk bottom) ---
bpy.ops.mesh.primitive_circle_add(vertices=32, radius=1, location=(0, 0, 0))

disk_obj = bpy.context.active_object
bpy.ops.object.editmode_toggle()
bpy.ops.mesh.fill()  # fills the circle to create a face
bpy.ops.object.editmode_toggle()

# --- 4. Join hemisphere and disk ---
bpy.ops.object.select_all(action='DESELECT')
hemi_obj.select_set(True)
disk_obj.select_set(True)
bpy.context.view_layer.objects.active = hemi_obj
bpy.ops.object.join()  # disk is now part of hemisphere mesh

# --- 5. Apply Solidify Modifier ---
solidify = hemi_obj.modifiers.new(name='Solidify', type='SOLIDIFY')
solidify.thickness = 0.01
solidify.offset = 0
solidify.use_even_offset = True

bpy.ops.object.modifier_apply(modifier=solidify.name)

print("✅ Hollow hemisphere with a flat bottom created successfully.")