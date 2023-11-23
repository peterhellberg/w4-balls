const std = @import("std");
const w4 = @import("w4");

const screen_width: f32 = @floatFromInt(w4.SCREEN_SIZE);
const screen_height: f32 = @floatFromInt(w4.SCREEN_SIZE);

// Coldfire GB Palette
// https://lospec.com/palette-list/coldfire-gb
const palette: [4]u32 = .{
    0x46425e, // Gray
    0x5b768d, // Green
    0xd17c7c, // Orange
    0xf6c6a8, // Yellow
};

const sqrt = std.math.sqrt;

var prng = std.rand.DefaultPrng.init(0);
const random = prng.random();

// 640 ought to be enough for anybody.
var memory: [640]u8 = undefined;
var fba = std.heap.FixedBufferAllocator.init(&memory);
const allocator = fba.allocator();

const ball_radius = 2;
const line_args = LineArgs{ .radius = 6 };

const Pair = @Vector(2, *Ball);

var lines = std.BoundedArray(Line, 8).init(0) catch {};
var balls = std.BoundedArray(Ball, 128).init(0) catch {};
var pairs = std.BoundedArray(Pair, 96).init(0) catch {};
var fakes = std.BoundedArray(Ball, 64).init(0) catch {};

var selected_ball: ?*Ball = null;
var selected_line: ?*Line = null;
var selected_line_start: bool = false;

var biggest: f32 = 0;

var mouse = Mouse{};
var gamepads = Gamepads{};

export fn start() void {
    init();
}

export fn update() void {
    input();
    simulate();
    draw();
}

fn fillCircle(x: f32, y: f32, r: f32, c: Color) void {
    w4.DRAW_COLORS.* = color(c);
    w4.oval(
        @intFromFloat(x - r),
        @intFromFloat(y - r),
        @intFromFloat(r * 2),
        @intFromFloat(r * 2),
    );
}

fn drawLine(x1: f32, y1: f32, x2: f32, y2: f32, c: Color) void {
    w4.DRAW_COLORS.* = color(c);
    w4.line(
        @intFromFloat(x1),
        @intFromFloat(y1),
        @intFromFloat(x2),
        @intFromFloat(y2),
    );
}

fn clear(c: u8) void {
    for (w4.FRAMEBUFFER) |*x| {
        x.* = c - 1 | (c - 1 << 2) | (c - 1 << 4) | (c - 1 << 6);
    }
}

fn pix(x: f32, y: f32, c: Color) void {
    w4.DRAW_COLORS.* = color(c);
    pixel(@intFromFloat(x), @intFromFloat(y));
}

fn pixel(x: i32, y: i32) void {
    if (x < 0 or x > 160 or y < 0 or y > 160) {
        return;
    }

    const ux: usize = @intCast(x);
    const uy: usize = @intCast(y);
    const idx: usize = (uy * 160 + ux) >> 2;
    const sx: u3 = @intCast(x);
    const shift = (sx & 0b11) * 2;
    const mask = @as(u8, 0b11) << shift;
    const palette_color: u8 = @intCast(w4.DRAW_COLORS.* & 0b1111);

    if (palette_color == 0) {
        return;
    }

    const c = (palette_color - 1) & 0b11;

    w4.FRAMEBUFFER[idx] = (c << shift) | (w4.FRAMEBUFFER[idx] & ~mask);
}

fn text(str: []const u8, x: i32, y: i32, c: Color) void {
    w4.DRAW_COLORS.* = color(c);
    w4.text(str, x, y);
}

fn init() void {
    w4.PALETTE.* = palette;

    for (0..100) |_| {
        const x = random.float(f32) * screen_width;
        const y = random.float(f32) * screen_height / 6;

        appendBall(x, y, .{ .radius = ball_radius });
    }

    appendBall(28, 33, .{ .radius = ball_radius * 3 });
    appendBall(28, 36, .{ .radius = ball_radius * 2 });

    addLine(80, 2, 90, 30, line_args); // Wiper
    addLine(-10, 20, 40, 100, line_args); // Left
    addLine(170, 20, 120, 100, line_args); // Right
    addLine(25, 120, 125, 140, line_args); // Floor
}

fn input() void {
    mouse.update();
    gamepads.update();

    if (mouse.rightReleased() or gamepads.b1(0) or gamepads.b2(0)) {
        var ball = newBall(mouse.x, mouse.y, .{
            .radius = 1 + (@round(random.float(f32) * 5)),
        });

        balls.append(ball) catch {};
    }

    if (mouse.leftPressed()) {
        selected_ball = null;
        for (balls.slice()) |*ball| {
            if (isPointInCircle(ball.px, ball.py, ball.radius, mouse.x, mouse.y)) {
                selected_ball = ball;
                break;
            }
        }

        selected_line = null;
        for (lines.slice()) |*line| {
            if (isPointInCircle(line.sx, line.sy, line.radius, mouse.x, mouse.y)) {
                selected_line = line;
                selected_line_start = true;
                break;
            }

            if (isPointInCircle(line.ex, line.ey, line.radius, mouse.x, mouse.y)) {
                selected_line = line;
                selected_line_start = false;
                break;
            }
        }
    }

    if (mouse.leftHeld()) {
        if (selected_line != null) {
            if (selected_line_start) {
                selected_line.?.sx = mouse.x;
                selected_line.?.sy = mouse.y;
            } else {
                selected_line.?.ex = mouse.x;
                selected_line.?.ey = mouse.y;
            }
        }
    }

    if (mouse.leftReleased()) {
        if (selected_ball != null) {
            selected_ball.?.vx = 5 * ((selected_ball.?.px) - mouse.x);
            selected_ball.?.vy = 5 * ((selected_ball.?.py) - mouse.y);
        }

        selected_ball = null;
        selected_line = null;
    }
}

fn simulate() void {
    const stable: f32 = 0.05;
    const sim_updates: usize = 4;
    const max_sim_steps: usize = 15;
    const sim_elapsed_time: f32 = 0.008;

    for (0..sim_updates) |_| {
        for (balls.slice()) |*ball| {
            ball.sim_time_remaining = sim_elapsed_time;
        }

        for (0..max_sim_steps) |_| {
            pairs.len = 0;
            fakes.len = 0;

            for (balls.slice()) |*ball| {
                if (ball.sim_time_remaining > 0.0) {
                    ball.ox = ball.px;
                    ball.oy = ball.py;

                    ball.ax = -ball.vx * 0.75;
                    ball.ay = -ball.vy * 0.75 + 100.0;

                    ball.vx += ball.ax * ball.sim_time_remaining;
                    ball.vy += ball.ay * ball.sim_time_remaining;

                    ball.px += ball.vx * ball.sim_time_remaining;
                    ball.py += ball.vy * ball.sim_time_remaining;

                    if (ball.px < 0) ball.px += screen_width;
                    if (ball.px >= screen_width) ball.px -= screen_width;
                    if (ball.py - ball.radius >= screen_height) ball.py -= screen_height + ball.radius;

                    if (ball.py < -100) ball.py = 80;
                    if (ball.py > 260) ball.py = 80;
                    if (ball.px < -100) ball.px = 80;
                    if (ball.px > 260) ball.px = 80;

                    if (@abs(ball.vx * ball.vx + ball.vy * ball.vy) < stable) {
                        ball.vx = 0;
                        ball.vy = 0;
                    }
                }
            }

            for (balls.slice()) |*ball| {
                for (lines.slice()) |*edge| {
                    const x1 = edge.ex - edge.sx;
                    const y1 = edge.ey - edge.sy;

                    const x2 = ball.px - edge.sx;
                    const y2 = ball.py - edge.sy;

                    const edge_length = x1 * x1 + y1 * y1;

                    const t = @max(0, @min(edge_length, (x1 * x2 + y1 * y2))) / edge_length;

                    const cx = edge.sx + t * x1;
                    const cy = edge.sy + t * y1;

                    const distance = sqrt(
                        (ball.px - cx) * (ball.px - cx) + (ball.py - cy) * (ball.py - cy),
                    );

                    if (distance <= (ball.radius + edge.radius) and !(ball.vx == 0 and ball.vy == 0)) {
                        var fakeball: *Ball = fakes.addOneAssumeCapacity();

                        fakeball.radius = edge.radius;
                        fakeball.mass = ball.mass * 0.8;
                        fakeball.px = cx;
                        fakeball.py = cy;
                        fakeball.vx = -ball.vx;
                        fakeball.vy = -ball.vy;

                        const overlap = 1.0 * (distance - ball.radius - fakeball.radius);

                        ball.px -= overlap * (ball.px - fakeball.px) / distance;
                        ball.py -= overlap * (ball.py - fakeball.py) / distance;

                        pairs.append(.{ ball, fakeball }) catch {};
                    }
                }

                for (balls.slice()) |*target| {
                    const overlaps = doCirclesOverlap(ball.px, ball.py, ball.radius, target.px, target.py, target.radius);

                    if (ball.px != target.px and overlaps and ball.vy != 0) {
                        const distance = sqrt(
                            (ball.px - target.px) * (ball.px - target.px) + (ball.py - target.py) * (ball.py - target.py),
                        );

                        const overlap = 0.5 * (distance - ball.radius - target.radius);

                        ball.px -= overlap * (ball.px - target.px) / distance;
                        ball.py -= overlap * (ball.py - target.py) / distance;

                        target.px += overlap * (ball.px - target.px) / distance;
                        target.py += overlap * (ball.py - target.py) / distance;

                        pairs.append(.{ ball, target }) catch {};
                    }
                }

                const intended_speed = sqrt(ball.vx * ball.vx + ball.vy * ball.vy);
                const actual_distance = sqrt(
                    (ball.px - ball.ox) * (ball.px - ball.ox) + (ball.py - ball.oy) * (ball.py - ball.oy),
                );
                const actual_time = actual_distance / intended_speed;

                ball.sim_time_remaining = ball.sim_time_remaining - actual_time;
            }

            const efficiency = 1.00;

            for (pairs.slice()) |pair| {
                var ball = pair[0];
                var target = pair[1];

                const distance = sqrt(
                    (ball.px - target.px) * (ball.px - target.px) + (ball.py - target.py) * (ball.py - target.py),
                );

                const nx = (target.px - ball.px) / distance;
                const ny = (target.py - ball.py) / distance;

                const tx = -ny;
                const ty = nx;

                const t1 = ball.vx * tx + ball.vy * ty;
                const t2 = target.vx * tx + target.vy * ty;

                const n1 = ball.vx * nx + ball.vy * ny;
                const n2 = target.vx * nx + target.vy * ny;

                const m1 = efficiency * (n1 * (ball.mass - target.mass) + 2.0 * target.mass * n2) / (ball.mass + target.mass);
                const m2 = efficiency * (n2 * (target.mass - ball.mass) + 2.0 * ball.mass * n1) / (ball.mass + target.mass);

                ball.vx = tx * t1 + nx * m1;
                ball.vy = ty * t1 + ny * m1;
                target.vx = tx * t2 + nx * m2;
                target.vy = ty * t2 + ny * m2;
            }

            for (pairs.slice()) |pair| {
                var ball = pair[0];
                var target = pair[1];

                if (@round(ball.radius) == @round(target.radius)) {
                    if (ball == target) continue;

                    ball.radius += 1;
                    ball.mass = ball.radius * 10;
                    ball.px = (ball.px + target.px) / 2;
                    ball.py = (ball.py + target.py) / 2;

                    target.radius = 0;
                    target.mass = 0;
                }
            }
        }
    }

    for (0..5) |_| {
        for (balls.slice(), 0..) |*ball, i| {
            if (ball.radius < 2) {
                _ = balls.orderedRemove(i);
                break;
            }
        }
    }

    biggest = 0;
    for (balls.slice()) |*ball| {
        if (ball.radius > biggest) biggest = ball.radius;
    }
}

fn draw() void {
    text("S", 1, 142, .Yellow);
    anyText(9, 142, .Green, .{@as(i32, @intFromFloat(biggest))});

    text("B", 1, 150, .Yellow);
    anyText(9, 150, .Orange, .{balls.len});

    for (lines.slice()) |*line| {
        var nx = -(line.ey - line.sy);
        var ny = (line.ex - line.sx);
        const d = sqrt((nx * nx) + (ny * ny));

        nx /= d;
        ny /= d;

        drawLine(
            line.sx + nx * line.radius * 0.8,
            line.sy + ny * line.radius * 0.8,
            line.ex + nx * line.radius * 0.8,
            line.ey + ny * line.radius * 0.8,
            .Green,
        );

        drawLine(
            line.sx,
            line.sy,
            line.ex,
            line.ey,
            .Green,
        );

        drawLine(
            line.sx - nx * line.radius * 0.8,
            line.sy - ny * line.radius * 0.8,
            line.ex - nx * line.radius * 0.8,
            line.ey - ny * line.radius * 0.8,
            .Green,
        );

        fillCircle(line.sx, line.sy, line.radius, .Green);
        fillCircle(line.sx, line.sy, line.radius - 1, .Green);
        fillCircle(line.ex, line.ey, line.radius, .Green);
        fillCircle(line.ex, line.ey, line.radius - 1, .Green);
    }

    for (balls.slice()) |*ball| {
        if (ball.radius > 1) {
            const bc: Color = ballColor(ball.radius);
            const ic: Color = ballColor(ball.radius - 1);
            const xc: Color = ballColor(ball.radius + 1);

            fillCircle(ball.px, ball.py, ball.radius, bc);
            fillCircle(ball.px, ball.py, ball.radius - 1, ic);
            fillCircle(ball.px, ball.py, ball.radius - 2, bc);
            fillCircle(ball.px, ball.py, 1, xc);
        }
    }

    if (selected_ball != null) {
        const ball = selected_ball.?;

        drawLine(ball.px, ball.py, mouse.x, mouse.y, .Green);
        fillCircle(ball.px, ball.py, ball.radius, .Green);
    }
}

fn newBall(x: f32, y: f32, args: BallArgs) Ball {
    return .{
        .px = x,
        .py = y,
        .radius = args.radius,
        .mass = args.radius * 10,
    };
}

fn appendBall(x: f32, y: f32, args: BallArgs) void {
    balls.append(newBall(x, y, args)) catch {};
}

fn addLine(sx: f32, sy: f32, ex: f32, ey: f32, args: LineArgs) void {
    lines.append(.{
        .sx = sx,
        .sy = sy,
        .ex = ex,
        .ey = ey,
        .radius = args.radius,
    }) catch {};
}

fn doCirclesOverlap(x1: f32, y1: f32, r1: f32, x2: f32, y2: f32, r2: f32) bool {
    return @abs(
        (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2),
    ) <= (r1 + r2) * (r1 + r2);
}

fn isPointInCircle(x1: f32, y1: f32, r1: f32, px: f32, py: f32) bool {
    return @abs(
        (x1 - px) * (x1 - px) + (y1 - py) * (y1 - py),
    ) < (r1 * r1);
}

fn color(c: Color) u3 {
    return @intFromEnum(c);
}

fn anyText(x: i32, y: i32, c: Color, args: anytype) void {
    const str = std.fmt.allocPrint(allocator, "{any}", args) catch "";
    defer allocator.free(str);

    text(str, x, y, c);
}

fn ballColor(f: f32) Color {
    return switch (@as(i32, @intFromFloat(f))) {
        2 => .Green,
        3 => .Orange,
        4 => .Yellow,
        5 => .Green,
        6 => .Orange,
        7 => .Yellow,
        8 => .Green,
        9 => .Orange,
        10 => .Yellow,
        else => .Green,
    };
}

const Color = enum {
    None,
    Gray,
    Green,
    Orange,
    Yellow,
};

const MouseData = struct {
    x: i16 = 0,
    y: i16 = 0,
    b: u8 = 0,
};

const Mouse = struct {
    x: f32 = 0,
    y: f32 = 0,

    data: MouseData = .{},
    prev: MouseData = .{},

    fn update(self: *Mouse) void {
        self.prev = self.data;

        self.data.x = w4.MOUSE_X.*;
        self.data.y = w4.MOUSE_Y.*;
        self.data.b = w4.MOUSE_BUTTONS.*;

        if (self.data.x >= 0 and self.data.x <= screen_width)
            self.x = @floatFromInt(@as(i32, @intCast(self.data.x)));

        if (self.data.y >= 0 and self.data.y <= screen_height)
            self.y = @floatFromInt(@as(i32, @intCast(self.data.y)));
    }

    fn leftPressed(self: *Mouse) bool {
        return (self.data.b & w4.MOUSE_LEFT != 0) and !(self.prev.b & w4.MOUSE_LEFT != 0);
    }

    fn leftHeld(self: *Mouse) bool {
        return (self.data.b & w4.MOUSE_LEFT != 0) and (self.prev.b & w4.MOUSE_LEFT != 0);
    }

    fn leftReleased(self: *Mouse) bool {
        return !(self.data.b & w4.MOUSE_LEFT != 0) and (self.prev.b & w4.MOUSE_LEFT != 0);
    }

    fn rightHeld(self: *Mouse) bool {
        return (self.data.b & w4.MOUSE_RIGHT != 0) and (self.prev.b & w4.MOUSE_RIGHT != 0);
    }

    fn rightReleased(self: *Mouse) bool {
        return !(self.data.b & w4.MOUSE_RIGHT != 0) and (self.prev.b & w4.MOUSE_RIGHT != 0);
    }
};

const Gamepads = struct {
    prev: [4]u8 = .{0} ** 4,
    data: [4]u8 = .{0} ** 4,

    fn update(self: *Gamepads) void {
        self.prev = self.data;

        self.data[0] = w4.GAMEPAD1.*;
        self.data[1] = w4.GAMEPAD2.*;
        self.data[2] = w4.GAMEPAD3.*;
        self.data[3] = w4.GAMEPAD4.*;
    }

    fn b1(self: *Gamepads, n: u2) bool {
        return self.released(n, w4.BUTTON_1);
    }

    fn b2(self: *Gamepads, n: u2) bool {
        return self.released(n, w4.BUTTON_2);
    }

    fn released(self: *Gamepads, n: u2, btn: u8) bool {
        return !(self.data[n] & btn != 0) and (self.prev[n] & btn != 0);
    }
};

const Ball = struct {
    px: f32 = 0,
    py: f32 = 0,

    vx: f32 = 0,
    vy: f32 = 0,

    ax: f32 = 0,
    ay: f32 = 0,

    ox: f32 = 0,
    oy: f32 = 0,

    radius: f32 = 0,
    mass: f32 = 0,
    friction: f32 = 0,

    sim_time_remaining: f32 = 0,
};

const BallArgs = struct {
    radius: f32 = 5.0,
};

const Line = struct {
    sx: f32 = 0,
    sy: f32 = 0,

    ex: f32 = 0,
    ey: f32 = 0,

    radius: f32 = 0,
};

const LineArgs = struct {
    radius: f32 = 4.0,
};
