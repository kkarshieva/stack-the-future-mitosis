import os
from flask import Flask, render_template, request, redirect, url_for, session, flash
from dotenv import load_dotenv
from supabase import create_client
from werkzeug.utils import secure_filename

load_dotenv()

#load keys
SUPABASE_URL = os.environ["SUPABASE_URL"]
SUPABASE_ANON_KEY = os.environ["SUPABASE_ANON_KEY"]
FLASK_SECRET_KEY = os.environ["FLASK_SECRET_KEY"]

#flask session
app = Flask(__name__)
app.secret_key = FLASK_SECRET_KEY

#helpers for supabase
def sb_public():
    return create_client(SUPABASE_URL, SUPABASE_ANON_KEY) 

def sb_user():
    sb = create_client(SUPABASE_URL, SUPABASE_ANON_KEY)
    sb.auth.set_session(session["access_token"], session["refresh_token"])
    return sb

def is_logged_in():
    return "access_token" in session and "refresh_token" in session

def require_login():
    if not is_logged_in():
        return redirect(url_for("login"))
    return None

def get_user_id():
    u = sb_user().auth.get_user()
    return u.user.id

UPLOAD_FOLDER = '/path/to/the/uploads'
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'}
app = Flask(__name__)
app.config['UPLOAD_FOLDER']=UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route("/", methods=["GET"])
def root():
    return render_template("base.html")


@app.route("/onboarding", methods=["GET", "POST"])
def onboarding():
    guard = require_login()
    if guard:
        return guard
    if request.method == "POST":
        request.data
    return render_template("onboarding.html")


@app.route("/matching", methods=["GET"])
def matching():
    guard = require_login()
    if guard:
        return guard
    return render_template("matching.html")


@app.route("/dashboard", methods=["GET"])
def dashboard():
    guard = require_login()
    if guard:
        return guard
    return render_template("dashboard.html")

#auth routes
@app.get("/signup")
def signup():
    return render_template("signup.html")

@app.post("/signup")
def signup_post(): #create user in supabase
    email = request.form.get("email", "").strip()
    password = request.form.get("password", "")

    if not email or not password:
        flash("Email and password required.")
        return redirect(url_for("signup"))

    try:
        resp = sb_public().auth.sign_up({"email": email, "password": password})
    except Exception as e:
        flash(f"Signup failed: {e}")
        return redirect(url_for("signup")) #error

    if getattr(resp, "session", None) is None:
        flash("Account created. Now log in (and confirm email if required).")
        return redirect(url_for("login")) #for authentictaed login later

    #store token in flask session cookies
    session["access_token"] = resp.session.access_token
    session["refresh_token"] = resp.session.refresh_token
    #redirect to onboarding
    return redirect(url_for("onboarding"))


@app.get("/login")
def login():
    return render_template("login.html")


@app.post("/login")
def login_post():
    email = request.form.get("email", "").strip()
    password = request.form.get("password", "")

    #authenticate
    if not email or not password:
        flash("Email and password required.")
        return redirect(url_for("login"))
    try:
        resp = sb_public().auth.sign_in_with_password({"email": email, "password": password})
    except Exception:
        flash("Login failed. Check email/password.")
        return redirect(url_for("login"))

    session["access_token"] = resp.session.access_token
    session["refresh_token"] = resp.session.refresh_token
    return redirect(url_for("onboarding"))


@app.post("/logout")
def logout():
    session.clear()
    return redirect(url_for("login"))
