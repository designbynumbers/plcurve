# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: plcTopology.proto

import sys
_b=sys.version_info[0]<3 and (lambda x:x) or (lambda x:x.encode('latin1'))
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
from google.protobuf import descriptor_pb2
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor.FileDescriptor(
  name='plcTopology.proto',
  package='',
  serialized_pb=_b('\n\x11plcTopology.proto\"\x86\x04\n\x06PDCode\x12\x0b\n\x03uid\x18\x01 \x01(\x04\x12\x0c\n\x04hash\x18\x06 \x01(\t\x12\x1b\n\x05\x65\x64ges\x18\x02 \x03(\x0b\x32\x0c.PDCode.Edge\x12#\n\tcrossings\x18\x03 \x03(\x0b\x32\x10.PDCode.Crossing\x12%\n\ncomponents\x18\x04 \x03(\x0b\x32\x11.PDCode.Component\x12\x1b\n\x05\x66\x61\x63\x65s\x18\x05 \x03(\x0b\x32\x0c.PDCode.Face\x1a\x44\n\x04\x45\x64ge\x12\x0c\n\x04head\x18\x01 \x01(\r\x12\x0f\n\x07headpos\x18\x02 \x01(\r\x12\x0c\n\x04tail\x18\x03 \x01(\r\x12\x0f\n\x07tailpos\x18\x04 \x01(\r\x1aG\n\x08\x43rossing\x12\x11\n\x05\x65\x64ges\x18\x01 \x03(\rB\x02\x10\x01\x12(\n\x04sign\x18\x02 \x01(\x0e\x32\x13.PDCode.Orientation:\x05UNSET\x1a+\n\tComponent\x12\x11\n\x05\x65\x64ges\x18\x01 \x03(\rB\x02\x10\x01\x12\x0b\n\x03tag\x18\x02 \x01(\r\x1ai\n\x04\x46\x61\x63\x65\x12$\n\x05\x65\x64ges\x18\x01 \x03(\x0b\x32\x15.PDCode.Face.FaceEdge\x1a;\n\x08\x46\x61\x63\x65\x45\x64ge\x12\x0c\n\x04\x65\x64ge\x18\x01 \x01(\r\x12!\n\x04sign\x18\x02 \x01(\x0e\x32\x13.PDCode.Orientation\"4\n\x0bOrientation\x12\x0c\n\x08NEGATIVE\x10\x00\x12\x0c\n\x08POSITIVE\x10\x01\x12\t\n\x05UNSET\x10\x02')
)
_sym_db.RegisterFileDescriptor(DESCRIPTOR)



_PDCODE_ORIENTATION = _descriptor.EnumDescriptor(
  name='Orientation',
  full_name='PDCode.Orientation',
  filename=None,
  file=DESCRIPTOR,
  values=[
    _descriptor.EnumValueDescriptor(
      name='NEGATIVE', index=0, number=0,
      options=None,
      type=None),
    _descriptor.EnumValueDescriptor(
      name='POSITIVE', index=1, number=1,
      options=None,
      type=None),
    _descriptor.EnumValueDescriptor(
      name='UNSET', index=2, number=2,
      options=None,
      type=None),
  ],
  containing_type=None,
  options=None,
  serialized_start=488,
  serialized_end=540,
)
_sym_db.RegisterEnumDescriptor(_PDCODE_ORIENTATION)


_PDCODE_EDGE = _descriptor.Descriptor(
  name='Edge',
  full_name='PDCode.Edge',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='head', full_name='PDCode.Edge.head', index=0,
      number=1, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='headpos', full_name='PDCode.Edge.headpos', index=1,
      number=2, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='tail', full_name='PDCode.Edge.tail', index=2,
      number=3, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='tailpos', full_name='PDCode.Edge.tailpos', index=3,
      number=4, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=193,
  serialized_end=261,
)

_PDCODE_CROSSING = _descriptor.Descriptor(
  name='Crossing',
  full_name='PDCode.Crossing',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='edges', full_name='PDCode.Crossing.edges', index=0,
      number=1, type=13, cpp_type=3, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=_descriptor._ParseOptions(descriptor_pb2.FieldOptions(), _b('\020\001'))),
    _descriptor.FieldDescriptor(
      name='sign', full_name='PDCode.Crossing.sign', index=1,
      number=2, type=14, cpp_type=8, label=1,
      has_default_value=True, default_value=2,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=263,
  serialized_end=334,
)

_PDCODE_COMPONENT = _descriptor.Descriptor(
  name='Component',
  full_name='PDCode.Component',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='edges', full_name='PDCode.Component.edges', index=0,
      number=1, type=13, cpp_type=3, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=_descriptor._ParseOptions(descriptor_pb2.FieldOptions(), _b('\020\001'))),
    _descriptor.FieldDescriptor(
      name='tag', full_name='PDCode.Component.tag', index=1,
      number=2, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=336,
  serialized_end=379,
)

_PDCODE_FACE_FACEEDGE = _descriptor.Descriptor(
  name='FaceEdge',
  full_name='PDCode.Face.FaceEdge',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='edge', full_name='PDCode.Face.FaceEdge.edge', index=0,
      number=1, type=13, cpp_type=3, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='sign', full_name='PDCode.Face.FaceEdge.sign', index=1,
      number=2, type=14, cpp_type=8, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=427,
  serialized_end=486,
)

_PDCODE_FACE = _descriptor.Descriptor(
  name='Face',
  full_name='PDCode.Face',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='edges', full_name='PDCode.Face.edges', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[_PDCODE_FACE_FACEEDGE, ],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=381,
  serialized_end=486,
)

_PDCODE = _descriptor.Descriptor(
  name='PDCode',
  full_name='PDCode',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='uid', full_name='PDCode.uid', index=0,
      number=1, type=4, cpp_type=4, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='hash', full_name='PDCode.hash', index=1,
      number=6, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='edges', full_name='PDCode.edges', index=2,
      number=2, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='crossings', full_name='PDCode.crossings', index=3,
      number=3, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='components', full_name='PDCode.components', index=4,
      number=4, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='faces', full_name='PDCode.faces', index=5,
      number=5, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[_PDCODE_EDGE, _PDCODE_CROSSING, _PDCODE_COMPONENT, _PDCODE_FACE, ],
  enum_types=[
    _PDCODE_ORIENTATION,
  ],
  options=None,
  is_extendable=False,
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=22,
  serialized_end=540,
)

_PDCODE_EDGE.containing_type = _PDCODE
_PDCODE_CROSSING.fields_by_name['sign'].enum_type = _PDCODE_ORIENTATION
_PDCODE_CROSSING.containing_type = _PDCODE
_PDCODE_COMPONENT.containing_type = _PDCODE
_PDCODE_FACE_FACEEDGE.fields_by_name['sign'].enum_type = _PDCODE_ORIENTATION
_PDCODE_FACE_FACEEDGE.containing_type = _PDCODE_FACE
_PDCODE_FACE.fields_by_name['edges'].message_type = _PDCODE_FACE_FACEEDGE
_PDCODE_FACE.containing_type = _PDCODE
_PDCODE.fields_by_name['edges'].message_type = _PDCODE_EDGE
_PDCODE.fields_by_name['crossings'].message_type = _PDCODE_CROSSING
_PDCODE.fields_by_name['components'].message_type = _PDCODE_COMPONENT
_PDCODE.fields_by_name['faces'].message_type = _PDCODE_FACE
_PDCODE_ORIENTATION.containing_type = _PDCODE
DESCRIPTOR.message_types_by_name['PDCode'] = _PDCODE

PDCode = _reflection.GeneratedProtocolMessageType('PDCode', (_message.Message,), dict(

  Edge = _reflection.GeneratedProtocolMessageType('Edge', (_message.Message,), dict(
    DESCRIPTOR = _PDCODE_EDGE,
    __module__ = 'plcTopology_pb2'
    # @@protoc_insertion_point(class_scope:PDCode.Edge)
    ))
  ,

  Crossing = _reflection.GeneratedProtocolMessageType('Crossing', (_message.Message,), dict(
    DESCRIPTOR = _PDCODE_CROSSING,
    __module__ = 'plcTopology_pb2'
    # @@protoc_insertion_point(class_scope:PDCode.Crossing)
    ))
  ,

  Component = _reflection.GeneratedProtocolMessageType('Component', (_message.Message,), dict(
    DESCRIPTOR = _PDCODE_COMPONENT,
    __module__ = 'plcTopology_pb2'
    # @@protoc_insertion_point(class_scope:PDCode.Component)
    ))
  ,

  Face = _reflection.GeneratedProtocolMessageType('Face', (_message.Message,), dict(

    FaceEdge = _reflection.GeneratedProtocolMessageType('FaceEdge', (_message.Message,), dict(
      DESCRIPTOR = _PDCODE_FACE_FACEEDGE,
      __module__ = 'plcTopology_pb2'
      # @@protoc_insertion_point(class_scope:PDCode.Face.FaceEdge)
      ))
    ,
    DESCRIPTOR = _PDCODE_FACE,
    __module__ = 'plcTopology_pb2'
    # @@protoc_insertion_point(class_scope:PDCode.Face)
    ))
  ,
  DESCRIPTOR = _PDCODE,
  __module__ = 'plcTopology_pb2'
  # @@protoc_insertion_point(class_scope:PDCode)
  ))
_sym_db.RegisterMessage(PDCode)
_sym_db.RegisterMessage(PDCode.Edge)
_sym_db.RegisterMessage(PDCode.Crossing)
_sym_db.RegisterMessage(PDCode.Component)
_sym_db.RegisterMessage(PDCode.Face)
_sym_db.RegisterMessage(PDCode.Face.FaceEdge)


_PDCODE_CROSSING.fields_by_name['edges'].has_options = True
_PDCODE_CROSSING.fields_by_name['edges']._options = _descriptor._ParseOptions(descriptor_pb2.FieldOptions(), _b('\020\001'))
_PDCODE_COMPONENT.fields_by_name['edges'].has_options = True
_PDCODE_COMPONENT.fields_by_name['edges']._options = _descriptor._ParseOptions(descriptor_pb2.FieldOptions(), _b('\020\001'))
# @@protoc_insertion_point(module_scope)